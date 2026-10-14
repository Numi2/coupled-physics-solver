#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "numi/matter/suture_thread.h"
#include <cmath>
#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>
int main(int argc,char** argv) { @autoreleasepool {try {
    if(argc!=2)throw std::runtime_error("usage: coupled_suture_contract METALLIB");
    auto device=MTLCreateSystemDefaultDevice();auto queue=[device newCommandQueue];NSError* error=nil;
    auto library=[device newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:argv[1]]] error:&error];
    if(!library)throw std::runtime_error("contract shader unavailable");
    auto pipeline=[device newComputePipelineStateWithFunction:[library newFunctionWithName:@"suture_ccd_contract"] error:&error];
    auto results=[device newBufferWithLength:9*sizeof(nm_float4) options:MTLResourceStorageModeShared];
    auto command=[queue commandBuffer];auto encoder=[command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];[encoder setBuffer:results offset:0 atIndex:0];
    [encoder dispatchThreads:MTLSizeMake(9,1,1) threadsPerThreadgroup:MTLSizeMake(9,1,1)];[encoder endEncoding];[command commit];[command waitUntilCompleted];
    const auto* r=static_cast<const nm_float4*>(results.contents);
    for(unsigned i=0;i<4;++i) {
        bool expectedHit=i==1 || i==3;
        if(r[i].y!=1 || r[i].x!=(expectedHit?1:0) || !std::isfinite(r[i].w))throw std::runtime_error("native CCD result violates analytic geometry");
        if(expectedHit && !(r[i].z>.49 && r[i].z<.50))throw std::runtime_error("CCD crossing bound outside expected interval");
    }
    for(unsigned i=4;i<7;++i) if(r[i].x!=1 || r[i].y!=0 || r[i].z!=0)throw std::runtime_error("small-scale mollifier erased a nondegenerate contact");
    if(std::abs(r[7].x-.00025f)>1e-9 || std::abs(r[7].y-.00025f)>1e-9 || std::abs(r[7].z)>1e-9 || std::abs(r[7].w-.25f)>1e-6)throw std::runtime_error("small triangle closest point is incorrect");
    if(std::abs(r[8].x-.5f)>1e-6 || std::abs(r[8].y-.5f)>1e-6 || std::abs(r[8].z-.00001f)>1e-9)throw std::runtime_error("small edge closest points are incorrect");
    // Compare the actual native thread tangent with central differences of
    // the actual native elastic force, using independently perturbed poses.
    std::array<NMFEMNodeStateGPU,4> nodes{};
    nodes[0].positionAndMass={0,-.0013f,.00075f,1};nodes[1].positionAndMass={0,.0013f,.00075f,1};
    nodes[2].positionAndMass={0,.0013f,-.00075f,1};nodes[3].positionAndMass={0,-.0013f,-.00075f,1};
    std::array<NMSutureAnchorGPU,4> anchors{};
    for(unsigned i=0;i<4;++i)anchors[i]={{i,0,i,0},{1,0,0,0}};
    NMSutureLoopGPU loop{{.0082f,.0142f,.5f,.0012f},{},{}};
    std::array<nm_float4,4> direction{{{.01f,.02f,-.01f,0},{-.02f,.01f,.015f,0},{.015f,-.01f,-.02f,0},{-.01f,-.02f,.01f,0}}};
    const auto make=[&](const void* bytes,NSUInteger length){return [device newBufferWithBytes:bytes length:length options:MTLResourceStorageModeShared];};
    auto nodeBuffer=make(nodes.data(),sizeof(nodes));auto anchorBuffer=make(anchors.data(),sizeof(anchors));auto loopBuffer=make(&loop,sizeof(loop));auto directionBuffer=make(direction.data(),sizeof(direction));
    auto forceBuffer=[device newBufferWithLength:sizeof(direction) options:MTLResourceStorageModeShared];
    std::array<nm_float4,4> zero{};auto tangentBuffer=make(zero.data(),sizeof(zero));
    auto forcePipeline=[device newComputePipelineStateWithFunction:[library newFunctionWithName:@"numi_matter_metal::nm_suture_thread_forces"] error:&error];
    auto operatorPipeline=[device newComputePipelineStateWithFunction:[library newFunctionWithName:@"numi_matter_metal::nm_suture_thread_operator"] error:&error];
    nm_uint4 counts{4,4,1,1};float progress=.5f,dt=.001f;
    const auto evaluate=[&](bool tangent){
        auto cb=[queue commandBuffer];auto e=[cb computeCommandEncoder];
        [e setComputePipelineState:forcePipeline];[e setBuffer:nodeBuffer offset:0 atIndex:0];[e setBuffer:anchorBuffer offset:0 atIndex:1];[e setBuffer:loopBuffer offset:0 atIndex:2];[e setBuffer:forceBuffer offset:0 atIndex:3];[e setBytes:&counts length:sizeof(counts) atIndex:4];[e setBytes:&progress length:sizeof(progress) atIndex:5];[e dispatchThreads:MTLSizeMake(4,1,1) threadsPerThreadgroup:MTLSizeMake(4,1,1)];
        if(tangent){[e setComputePipelineState:operatorPipeline];[e setBuffer:directionBuffer offset:0 atIndex:3];[e setBuffer:tangentBuffer offset:0 atIndex:4];[e setBytes:&counts length:sizeof(counts) atIndex:5];[e setBytes:&dt length:sizeof(dt) atIndex:6];[e dispatchThreads:MTLSizeMake(4,1,1) threadsPerThreadgroup:MTLSizeMake(4,1,1)];}
        [e endEncoding];[cb commit];[cb waitUntilCompleted];
        if(cb.status!=MTLCommandBufferStatusCompleted)throw std::runtime_error("native thread contract command failed");
    };
    evaluate(true);std::array<nm_float4,4> tangent{};std::memcpy(tangent.data(),tangentBuffer.contents,sizeof(tangent));
    std::array<std::array<nm_float4,4>,2> perturbed{};const float epsilon=1e-4f;
    for(unsigned sign=0;sign<2;++sign){auto moved=nodes;
        for(unsigned i=0;i<4;++i){auto& p=moved[i].positionAndMass;const auto d=direction[i];float h=(sign?1:-1)*epsilon;p.x+=h*d.x;p.y+=h*d.y;p.z+=h*d.z;}
        std::memcpy(nodeBuffer.contents,moved.data(),sizeof(moved));evaluate(false);std::memcpy(perturbed[sign].data(),forceBuffer.contents,sizeof(direction));
    }
    double errorSquared=0,referenceSquared=0;
    for(unsigned i=0;i<4;++i)for(unsigned axis=0;axis<3;++axis){
        const float* a=reinterpret_cast<const float*>(&perturbed[0][i]);const float* b=reinterpret_cast<const float*>(&perturbed[1][i]);const float* t=reinterpret_cast<const float*>(&tangent[i]);
        double expected=-dt*dt*(b[axis]-a[axis])/(2*epsilon);errorSquared+=std::pow(t[axis]-expected,2);referenceSquared+=expected*expected;
    }
    double tangentError=std::sqrt(errorSquared/std::max(referenceSquared,1e-30));
    if(!(tangentError<.002))throw std::runtime_error("native thread tangent disagrees with force finite differences");
    std::cout<<"suture_thread_tangent=ok relative_error="<<tangentError<<'\n';
    std::cout<<"suture_ccd_contract=ok cases=9 closest_points=true scale_invariance=true translation_invariance=true crossing_bounds=true\n";
    return 0;
} catch(const std::exception& e) {std::cerr<<"suture_ccd_contract=failed "<<e.what()<<'\n';return 1;}} }
