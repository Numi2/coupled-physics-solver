# Matter, contact, and suture research

Coupled Physics Solver owns the Matter, coupled stepping, FEM, contact, and suture research in this handoff. The ledger includes the former Numi Lab `coupled` integration branch. A compatibility tag named `coupled` preserves its final revision for existing pinned checkout instructions.

This handoff indexes 52 former Numi Lab branch heads and imports 204 byte-preserved source and evidence files. The snapshots are historical research inputs stored under `docs` or `Docs`. Reproducing a result requires the original source revision, hardware, inputs, and acceptance criteria recorded in its history. The owning repository’s live build keeps its existing source and targets.

## Branch and file records

- [Branch ledger](branches.tsv) records the original branch name, tip, root tree, archive tag, comparison base, and selected snapshot paths.
- [Changed paths](changed-paths.tsv.gz) lists every path changed against the ledger’s comparison base, including deletions and original Git object IDs. Branches that share their tip with their comparison base have an empty delta, and branches without a shared base are compared with an empty tree.
- [Source manifest](source-files.tsv) and [SHA-256 checksums](source-files.sha256) record the original revision, source path, Git mode, blob ID, and content of each imported file. From this directory, run `shasum -a 256 -c source-files.sha256` to verify the imports.

## Original Numi Lab draft reviews

- [Matter compiler foundation (Numi Lab #2)](https://github.com/Numi2/numi-lab/pull/2). The draft head is an ancestor of Numi Lab `main` and the final `coupled` revision. This package preserves the full review diff against the draft’s original target and records the original tip in the branch ledger. [Original diff](source-reviews/numi-lab-pr-2.diff).

[Review diff checksums](source-reviews/SHA256SUMS).

## Restore a branch

The complete Git history and bulk artifacts remain available from the annotated archive tags in [Numi Lab](https://github.com/Numi2/numi-lab/tags). Copy the original branch name and archive tag from the ledger. This example fetches the complete object history for that tag into a fresh local repository and recreates the original branch.

```sh
git init --bare restored-numi-lab.git
branch='<source_branch from branches.tsv>'
tag='<archive_tag from branches.tsv>'
git -C restored-numi-lab.git fetch --no-tags https://github.com/Numi2/numi-lab.git \
  "refs/tags/$tag:refs/tags/$tag"
git -C restored-numi-lab.git branch "$branch" "refs/tags/$tag^{commit}"
git -C restored-numi-lab.git rev-parse "refs/heads/$branch" "refs/heads/$branch^{tree}"
```

The final command prints the source tip and root tree from the ledger.
