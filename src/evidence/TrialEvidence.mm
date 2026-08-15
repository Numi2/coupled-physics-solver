#import <Foundation/Foundation.h>

#include "coupled/TrialEvidence.hpp"

#include <CommonCrypto/CommonDigest.h>

#include <array>
#include <bit>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <system_error>

#include <unistd.h>

namespace coupled {
namespace {

constexpr std::uint64_t kFnvOffset = 1469598103934665603ull;
constexpr std::uint64_t kFnvPrime = 1099511628211ull;

EvidenceDiagnostics accept() {
    return {.accepted = true, .message = "accepted"};
}

EvidenceDiagnostics reject(std::string message) {
    return {.accepted = false, .message = std::move(message)};
}

bool finite(const double value) noexcept {
    return std::isfinite(value);
}

bool validSha256(const std::string_view value) noexcept {
    if (value.size() != 64u) {
        return false;
    }
    for (const char character : value) {
        if (!((character >= '0' && character <= '9') ||
              (character >= 'a' && character <= 'f'))) {
            return false;
        }
    }
    return true;
}

bool validRevision(const std::string_view value) noexcept {
    if (value.size() != 40u) {
        return false;
    }
    for (const char character : value) {
        if (!((character >= '0' && character <= '9') ||
              (character >= 'a' && character <= 'f'))) {
            return false;
        }
    }
    return true;
}

bool safeRelativePath(const std::string& value) {
    if (value.empty()) {
        return false;
    }
    const std::filesystem::path path(value);
    if (path.is_absolute() || path.has_root_path()) {
        return false;
    }
    const auto normalized = path.lexically_normal();
    if (normalized.empty() || normalized == ".") {
        return false;
    }
    for (const auto& component : normalized) {
        if (component == "..") {
            return false;
        }
    }
    return normalized.generic_string() == path.generic_string();
}

std::uint32_t scalarWidth(const std::string_view value) noexcept {
    if (value == "f32le" || value == "u32le") {
        return 4u;
    }
    if (value == "f64le" || value == "i64le") {
        return 8u;
    }
    if (value == "u8") {
        return 1u;
    }
    return 0u;
}

std::string escaped(const std::string_view value) {
    std::ostringstream output;
    output << '"';
    for (const unsigned char character : value) {
        switch (character) {
        case '"': output << "\\\""; break;
        case '\\': output << "\\\\"; break;
        case '\b': output << "\\b"; break;
        case '\f': output << "\\f"; break;
        case '\n': output << "\\n"; break;
        case '\r': output << "\\r"; break;
        case '\t': output << "\\t"; break;
        default:
            if (character < 0x20u) {
                output << "\\u" << std::hex << std::setw(4)
                       << std::setfill('0')
                       << static_cast<unsigned int>(character)
                       << std::dec;
            } else {
                output << static_cast<char>(character);
            }
        }
    }
    output << '"';
    return output.str();
}

std::string hex64(const std::uint64_t value) {
    std::ostringstream output;
    output << std::hex << std::setw(16) << std::setfill('0') << value;
    return output.str();
}

bool parseHex64(NSString* value, std::uint64_t& output) {
    if (![value isKindOfClass:[NSString class]] || value.length != 16u) {
        return false;
    }
    const std::string text(value.UTF8String == nullptr ? "" : value.UTF8String);
    std::uint64_t parsed = 0u;
    for (const char character : text) {
        std::uint8_t digit = 0u;
        if (character >= '0' && character <= '9') {
            digit = static_cast<std::uint8_t>(character - '0');
        } else if (character >= 'a' && character <= 'f') {
            digit = static_cast<std::uint8_t>(character - 'a' + 10);
        } else {
            return false;
        }
        parsed = (parsed << 4u) | digit;
    }
    output = parsed;
    return true;
}

void writeContentReference(
    std::ostream& output,
    const ContentReference& reference
) {
    output << "{\"byteCount\":" << reference.byteCount
           << ",\"path\":" << escaped(reference.path)
           << ",\"sha256\":" << escaped(reference.sha256) << '}';
}

NSString* stringField(NSDictionary* dictionary, NSString* key) {
    id value = dictionary[key];
    return [value isKindOfClass:[NSString class]] ? value : nil;
}

NSNumber* numberField(NSDictionary* dictionary, NSString* key) {
    id value = dictionary[key];
    return [value isKindOfClass:[NSNumber class]] ? value : nil;
}

bool readContentReference(id value, ContentReference& output) {
    if (![value isKindOfClass:[NSDictionary class]]) {
        return false;
    }
    NSDictionary* dictionary = value;
    NSString* path = stringField(dictionary, @"path");
    NSString* hash = stringField(dictionary, @"sha256");
    NSNumber* bytes = numberField(dictionary, @"byteCount");
    if (path == nil || hash == nil || bytes == nil || bytes.longLongValue < 0) {
        return false;
    }
    output.path = path.UTF8String == nullptr ? "" : path.UTF8String;
    output.sha256 = hash.UTF8String == nullptr ? "" : hash.UTF8String;
    output.byteCount = bytes.unsignedLongLongValue;
    return true;
}

std::optional<CohortSplit> parseCohortSplit(NSString* value) {
    if (value == nil) {
        return std::nullopt;
    }
    const std::string_view text(value.UTF8String == nullptr ? "" : value.UTF8String);
    if (text == "calibration") return CohortSplit::calibration;
    if (text == "model-selection") return CohortSplit::modelSelection;
    if (text == "blind-transfer") return CohortSplit::blindTransfer;
    return std::nullopt;
}

std::optional<SpecimenRole> parseSpecimenRole(NSString* value) {
    if (value == nil) {
        return std::nullopt;
    }
    const std::string_view text(value.UTF8String == nullptr ? "" : value.UTF8String);
    if (text == "bench-identification") return SpecimenRole::benchIdentification;
    if (text == "dvrk-closure") return SpecimenRole::dvrkClosure;
    return std::nullopt;
}

std::uint64_t hashBytes(
    std::uint64_t hash,
    const void* bytes,
    const std::size_t count
) noexcept {
    const auto* values = static_cast<const std::uint8_t*>(bytes);
    for (std::size_t index = 0u; index < count; ++index) {
        hash ^= values[index];
        hash *= kFnvPrime;
    }
    return hash;
}

} // namespace

std::string_view cohortSplitName(const CohortSplit value) noexcept {
    switch (value) {
    case CohortSplit::calibration: return "calibration";
    case CohortSplit::modelSelection: return "model-selection";
    case CohortSplit::blindTransfer: return "blind-transfer";
    }
    return "invalid";
}

std::string_view specimenRoleName(const SpecimenRole value) noexcept {
    switch (value) {
    case SpecimenRole::benchIdentification: return "bench-identification";
    case SpecimenRole::dvrkClosure: return "dvrk-closure";
    }
    return "invalid";
}

EvidenceDiagnostics validateTrialManifest(
    const CoupledTrialManifest& manifest
) {
    if (manifest.schemaVersion != kTrialManifestSchemaVersion) {
        return reject("unsupported trial manifest schema version");
    }
    if (manifest.trialId.empty() || manifest.animalId.empty() ||
        manifest.specimenId.empty() || manifest.pairedSpecimenId.empty()) {
        return reject("trial, animal, specimen, and paired-specimen identifiers are required");
    }
    if (manifest.specimenId == manifest.pairedSpecimenId) {
        return reject("paired specimens must have distinct identifiers");
    }
    if (cohortSplitName(manifest.cohortSplit) == "invalid" ||
        specimenRoleName(manifest.specimenRole) == "invalid") {
        return reject("trial cohort split or specimen role is invalid");
    }
    if (manifest.captureStartUnixNs <= 0) {
        return reject("capture start must use a positive Unix nanosecond timestamp");
    }
    if (!validRevision(manifest.sourceRevision)) {
        return reject("source revision must be a lowercase 40-character Git hash");
    }
    if (manifest.worldFingerprint == 0u ||
        manifest.deviceProgramFingerprint == 0u) {
        return reject("world and device-program fingerprints are required");
    }
    if (manifest.calibrations.empty()) {
        return reject("at least one rig calibration is required");
    }
    if (manifest.streams.empty()) {
        return reject("at least one trial stream is required");
    }

    std::set<std::string> paths;
    for (const auto& calibration : manifest.calibrations) {
        if (!safeRelativePath(calibration.path) ||
            !validSha256(calibration.sha256) ||
            calibration.byteCount == 0u) {
            return reject("calibration references require a normalized relative path, SHA-256, and nonzero size");
        }
        if (!paths.insert(calibration.path).second) {
            return reject("trial manifest contains a duplicate payload path");
        }
    }

    std::set<std::string> names;
    for (const auto& stream : manifest.streams) {
        const std::uint32_t width = scalarWidth(stream.scalarType);
        if (stream.name.empty() || stream.clockId.empty() || width == 0u ||
            stream.channelCount == 0u || stream.sampleCount == 0u ||
            !finite(stream.nominalRateHz) || stream.nominalRateHz <= 0.0) {
            return reject("trial streams require a unique name, clock, supported scalar type, dimensions, and finite positive rate");
        }
        const std::uint64_t packedWidth =
            static_cast<std::uint64_t>(width) * stream.channelCount;
        if (stream.strideBytes < packedWidth ||
            stream.content.byteCount !=
                stream.sampleCount * static_cast<std::uint64_t>(stream.strideBytes)) {
            return reject("trial stream stride or byte count is inconsistent with its shape");
        }
        if (!safeRelativePath(stream.content.path) ||
            !validSha256(stream.content.sha256)) {
            return reject("trial stream payload reference is invalid");
        }
        if (!names.insert(stream.name).second) {
            return reject("trial stream names must be unique");
        }
        if (!paths.insert(stream.content.path).second) {
            return reject("trial manifest contains a duplicate payload path");
        }
    }
    return accept();
}

std::string canonicalTrialManifestJSON(
    const CoupledTrialManifest& manifest
) {
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::setprecision(17);
    output << "{\"animalId\":" << escaped(manifest.animalId)
           << ",\"calibrations\":[";
    for (std::size_t index = 0u; index < manifest.calibrations.size(); ++index) {
        if (index != 0u) output << ',';
        writeContentReference(output, manifest.calibrations[index]);
    }
    output << "],\"captureStartUnixNs\":" << manifest.captureStartUnixNs
           << ",\"cohortSplit\":" << escaped(cohortSplitName(manifest.cohortSplit))
           << ",\"deviceProgramFingerprint\":"
           << escaped(hex64(manifest.deviceProgramFingerprint))
           << ",\"pairedSpecimenId\":" << escaped(manifest.pairedSpecimenId)
           << ",\"schemaVersion\":" << manifest.schemaVersion
           << ",\"sourceRevision\":" << escaped(manifest.sourceRevision)
           << ",\"specimenId\":" << escaped(manifest.specimenId)
           << ",\"specimenRole\":" << escaped(specimenRoleName(manifest.specimenRole))
           << ",\"streams\":[";
    for (std::size_t index = 0u; index < manifest.streams.size(); ++index) {
        if (index != 0u) output << ',';
        const auto& stream = manifest.streams[index];
        output << "{\"channelCount\":" << stream.channelCount
               << ",\"clockId\":" << escaped(stream.clockId)
               << ",\"content\":";
        writeContentReference(output, stream.content);
        output << ",\"name\":" << escaped(stream.name)
               << ",\"nominalRateHz\":" << stream.nominalRateHz
               << ",\"sampleCount\":" << stream.sampleCount
               << ",\"scalarType\":" << escaped(stream.scalarType)
               << ",\"startTimeNs\":" << stream.startTimeNs
               << ",\"strideBytes\":" << stream.strideBytes << '}';
    }
    output << "],\"trialId\":" << escaped(manifest.trialId)
           << ",\"worldFingerprint\":" << escaped(hex64(manifest.worldFingerprint))
           << "}\n";
    return output.str();
}

EvidenceDiagnostics writeTrialManifest(
    const std::filesystem::path& path,
    const CoupledTrialManifest& manifest
) {
    if (const auto diagnostics = validateTrialManifest(manifest);
        !diagnostics.accepted) {
        return diagnostics;
    }
    std::error_code filesystemError;
    if (std::filesystem::exists(path, filesystemError) || filesystemError) {
        return reject("trial manifests are immutable and the destination already exists");
    }
    const auto temporary = path.string() + ".tmp." +
        std::to_string(static_cast<unsigned long long>(::getpid()));
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
        return reject("failed to open trial manifest for writing");
    }
    const std::string encoded = canonicalTrialManifestJSON(manifest);
    output.write(encoded.data(), static_cast<std::streamsize>(encoded.size()));
    output.close();
    if (!output) {
        std::filesystem::remove(temporary, filesystemError);
        return reject("failed to write complete trial manifest");
    }
    std::filesystem::rename(temporary, path, filesystemError);
    if (filesystemError) {
        std::filesystem::remove(temporary, filesystemError);
        return reject("failed to publish trial manifest atomically");
    }
    return accept();
}

EvidenceDiagnostics readTrialManifest(
    const std::filesystem::path& path,
    CoupledTrialManifest& manifest
) {
    @autoreleasepool {
        NSData* data = [NSData dataWithContentsOfFile:
            [NSString stringWithUTF8String:path.c_str()]];
        if (data == nil) {
            return reject("failed to read trial manifest");
        }
        NSError* error = nil;
        id root = [NSJSONSerialization JSONObjectWithData:data options:0 error:&error];
        if (root == nil || ![root isKindOfClass:[NSDictionary class]]) {
            const char* message = error == nil ? nullptr : error.localizedDescription.UTF8String;
            return reject(message == nullptr ? "invalid trial manifest JSON" : message);
        }
        NSDictionary* dictionary = root;
        CoupledTrialManifest staged;
        NSNumber* schemaVersion = numberField(dictionary, @"schemaVersion");
        NSNumber* captureStart = numberField(dictionary, @"captureStartUnixNs");
        NSString* trialId = stringField(dictionary, @"trialId");
        NSString* animalId = stringField(dictionary, @"animalId");
        NSString* specimenId = stringField(dictionary, @"specimenId");
        NSString* pairedId = stringField(dictionary, @"pairedSpecimenId");
        NSString* sourceRevision = stringField(dictionary, @"sourceRevision");
        const auto split = parseCohortSplit(stringField(dictionary, @"cohortSplit"));
        const auto role = parseSpecimenRole(stringField(dictionary, @"specimenRole"));
        if (schemaVersion == nil || captureStart == nil || trialId == nil ||
            animalId == nil || specimenId == nil || pairedId == nil ||
            sourceRevision == nil || !split || !role ||
            !parseHex64(stringField(dictionary, @"worldFingerprint"), staged.worldFingerprint) ||
            !parseHex64(stringField(dictionary, @"deviceProgramFingerprint"),
                        staged.deviceProgramFingerprint)) {
            return reject("trial manifest has missing or invalid scalar fields");
        }
        staged.schemaVersion = schemaVersion.unsignedIntValue;
        staged.captureStartUnixNs = captureStart.longLongValue;
        staged.trialId = trialId.UTF8String == nullptr ? "" : trialId.UTF8String;
        staged.animalId = animalId.UTF8String == nullptr ? "" : animalId.UTF8String;
        staged.specimenId = specimenId.UTF8String == nullptr ? "" : specimenId.UTF8String;
        staged.pairedSpecimenId = pairedId.UTF8String == nullptr ? "" : pairedId.UTF8String;
        staged.sourceRevision = sourceRevision.UTF8String == nullptr ? "" : sourceRevision.UTF8String;
        staged.cohortSplit = *split;
        staged.specimenRole = *role;

        id calibrationValues = dictionary[@"calibrations"];
        id streamValues = dictionary[@"streams"];
        if (![calibrationValues isKindOfClass:[NSArray class]] ||
            ![streamValues isKindOfClass:[NSArray class]]) {
            return reject("trial manifest calibration and stream fields must be arrays");
        }
        for (id value in static_cast<NSArray*>(calibrationValues)) {
            ContentReference reference;
            if (!readContentReference(value, reference)) {
                return reject("trial manifest contains an invalid calibration reference");
            }
            staged.calibrations.push_back(std::move(reference));
        }
        for (id value in static_cast<NSArray*>(streamValues)) {
            if (![value isKindOfClass:[NSDictionary class]]) {
                return reject("trial manifest contains a non-object stream");
            }
            NSDictionary* streamDictionary = value;
            TrialStream stream;
            NSString* name = stringField(streamDictionary, @"name");
            NSString* type = stringField(streamDictionary, @"scalarType");
            NSString* clock = stringField(streamDictionary, @"clockId");
            NSNumber* channels = numberField(streamDictionary, @"channelCount");
            NSNumber* stride = numberField(streamDictionary, @"strideBytes");
            NSNumber* count = numberField(streamDictionary, @"sampleCount");
            NSNumber* start = numberField(streamDictionary, @"startTimeNs");
            NSNumber* rate = numberField(streamDictionary, @"nominalRateHz");
            if (name == nil || type == nil || clock == nil || channels == nil ||
                stride == nil || count == nil || start == nil || rate == nil ||
                !readContentReference(streamDictionary[@"content"], stream.content)) {
                return reject("trial manifest contains an invalid stream record");
            }
            stream.name = name.UTF8String == nullptr ? "" : name.UTF8String;
            stream.scalarType = type.UTF8String == nullptr ? "" : type.UTF8String;
            stream.clockId = clock.UTF8String == nullptr ? "" : clock.UTF8String;
            stream.channelCount = channels.unsignedIntValue;
            stream.strideBytes = stride.unsignedIntValue;
            stream.sampleCount = count.unsignedLongLongValue;
            stream.startTimeNs = start.longLongValue;
            stream.nominalRateHz = rate.doubleValue;
            staged.streams.push_back(std::move(stream));
        }
        if (const auto diagnostics = validateTrialManifest(staged);
            !diagnostics.accepted) {
            return diagnostics;
        }
        manifest = std::move(staged);
        return accept();
    }
}

std::string sha256File(
    const std::filesystem::path& path,
    std::string* error
) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        if (error != nullptr) *error = "failed to open payload for hashing";
        return {};
    }
    CC_SHA256_CTX context;
    CC_SHA256_Init(&context);
    std::array<char, 1u << 16u> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto count = input.gcount();
        if (count > 0) {
            CC_SHA256_Update(
                &context,
                buffer.data(),
                static_cast<CC_LONG>(count)
            );
        }
    }
    if (!input.eof()) {
        if (error != nullptr) *error = "failed while reading payload for hashing";
        return {};
    }
    std::array<unsigned char, CC_SHA256_DIGEST_LENGTH> digest{};
    CC_SHA256_Final(digest.data(), &context);
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (const unsigned char value : digest) {
        output << std::setw(2) << static_cast<unsigned int>(value);
    }
    if (error != nullptr) error->clear();
    return output.str();
}

EvidenceDiagnostics verifyTrialPayloads(
    const std::filesystem::path& manifestDirectory,
    const CoupledTrialManifest& manifest
) {
    if (const auto diagnostics = validateTrialManifest(manifest);
        !diagnostics.accepted) {
        return diagnostics;
    }
    const auto verify = [&](const ContentReference& reference) {
        const auto path = manifestDirectory / reference.path;
        std::error_code filesystemError;
        const auto byteCount = std::filesystem::file_size(path, filesystemError);
        if (filesystemError || byteCount != reference.byteCount) {
            return reject("trial payload size mismatch: " + reference.path);
        }
        std::string hashError;
        const std::string hash = sha256File(path, &hashError);
        if (!hashError.empty() || hash != reference.sha256) {
            return reject("trial payload SHA-256 mismatch: " + reference.path);
        }
        return accept();
    };
    for (const auto& calibration : manifest.calibrations) {
        if (const auto diagnostics = verify(calibration); !diagnostics.accepted) {
            return diagnostics;
        }
    }
    for (const auto& stream : manifest.streams) {
        if (const auto diagnostics = verify(stream.content); !diagnostics.accepted) {
            return diagnostics;
        }
    }
    return accept();
}

std::uint64_t trialManifestFingerprint(
    const CoupledTrialManifest& manifest
) {
    const std::string canonical = canonicalTrialManifestJSON(manifest);
    const std::uint64_t hash = hashBytes(
        kFnvOffset,
        canonical.data(),
        canonical.size()
    );
    return hash == 0u ? 1u : hash;
}

EvidenceDiagnostics validatePhysicsEvidence(
    const PhysicsEvidence& evidence
) {
    if (evidence.schemaVersion != kPhysicsEvidenceSchemaVersion) {
        return reject("unsupported physics evidence schema version");
    }
    if (evidence.evidenceId.empty() || evidence.deviceName.empty() ||
        evidence.operatingSystem.empty() || !validRevision(evidence.sourceRevision) ||
        !validSha256(evidence.trialManifestSha256) ||
        !validSha256(evidence.metallibSha256)) {
        return reject("physics evidence provenance is incomplete");
    }
    if (evidence.worldFingerprint == 0u ||
        evidence.deviceProgramFingerprint == 0u ||
        evidence.replayFingerprint == 0u) {
        return reject("physics evidence fingerprints are incomplete");
    }
    if (!finite(evidence.elapsedSeconds) || evidence.elapsedSeconds <= 0.0 ||
        !finite(evidence.environmentStepsPerSecond) ||
        evidence.environmentStepsPerSecond < 0.0 ||
        !finite(evidence.physicsSubstepsPerSecond) ||
        evidence.physicsSubstepsPerSecond < 0.0) {
        return reject("physics evidence timing is invalid");
    }
    if (!evidence.exactReplay || !evidence.fp64Parity ||
        !evidence.physicalOutcomesAccepted || evidence.failedSteps != 0u) {
        return reject("an absolute physics promotion gate failed");
    }
    if (!evidence.zeroSwapGrowth || !evidence.timelineCaptured) {
        return reject("memory or Metal timeline promotion evidence is missing");
    }
    if (!evidence.gpuCountersCaptured && evidence.gpuCounterStatus.empty()) {
        return reject("unsupported GPU counters must be recorded explicitly");
    }
    if (evidence.physicalMetrics.empty()) {
        return reject("physics evidence requires physical outcome metrics");
    }
    for (const auto& metric : evidence.physicalMetrics) {
        if (metric.name.empty() || metric.unit.empty() ||
            !finite(metric.observed) || !finite(metric.reference) ||
            !finite(metric.absoluteError) || !finite(metric.acceptanceBound) ||
            metric.absoluteError < 0.0 || metric.acceptanceBound < 0.0 ||
            !metric.accepted || metric.absoluteError > metric.acceptanceBound) {
            return reject("a physical outcome metric is invalid or outside its bound");
        }
    }
    return accept();
}

std::string canonicalPhysicsEvidenceJSON(
    const PhysicsEvidence& evidence
) {
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::boolalpha << std::setprecision(17)
           << "{\"deviceName\":" << escaped(evidence.deviceName)
           << ",\"deviceProgramFingerprint\":"
           << escaped(hex64(evidence.deviceProgramFingerprint))
           << ",\"elapsedSeconds\":" << evidence.elapsedSeconds
           << ",\"environmentSteps\":" << evidence.environmentSteps
           << ",\"environmentStepsPerSecond\":" << evidence.environmentStepsPerSecond
           << ",\"evidenceId\":" << escaped(evidence.evidenceId)
           << ",\"exactReplay\":" << evidence.exactReplay
           << ",\"failedSteps\":" << evidence.failedSteps
           << ",\"fp64Parity\":" << evidence.fp64Parity
           << ",\"gpuCounterStatus\":" << escaped(evidence.gpuCounterStatus)
           << ",\"gpuCountersCaptured\":" << evidence.gpuCountersCaptured
           << ",\"metallibSha256\":" << escaped(evidence.metallibSha256)
           << ",\"operatingSystem\":" << escaped(evidence.operatingSystem)
           << ",\"peakBytes\":" << evidence.peakBytes
           << ",\"physicalMetrics\":[";
    for (std::size_t index = 0u; index < evidence.physicalMetrics.size(); ++index) {
        if (index != 0u) output << ',';
        const auto& metric = evidence.physicalMetrics[index];
        output << "{\"absoluteError\":" << metric.absoluteError
               << ",\"acceptanceBound\":" << metric.acceptanceBound
               << ",\"accepted\":" << metric.accepted
               << ",\"name\":" << escaped(metric.name)
               << ",\"observed\":" << metric.observed
               << ",\"reference\":" << metric.reference
               << ",\"unit\":" << escaped(metric.unit) << '}';
    }
    output << "],\"physicalOutcomesAccepted\":"
           << evidence.physicalOutcomesAccepted
           << ",\"physicsSubsteps\":" << evidence.physicsSubsteps
           << ",\"physicsSubstepsPerSecond\":" << evidence.physicsSubstepsPerSecond
           << ",\"replayFingerprint\":" << escaped(hex64(evidence.replayFingerprint))
           << ",\"retainedBytes\":" << evidence.retainedBytes
           << ",\"schemaVersion\":" << evidence.schemaVersion
           << ",\"sourceRevision\":" << escaped(evidence.sourceRevision)
           << ",\"timelineCaptured\":" << evidence.timelineCaptured
           << ",\"trialManifestSha256\":" << escaped(evidence.trialManifestSha256)
           << ",\"worldFingerprint\":" << escaped(hex64(evidence.worldFingerprint))
           << ",\"zeroSwapGrowth\":" << evidence.zeroSwapGrowth << "}\n";
    return output.str();
}

} // namespace coupled
