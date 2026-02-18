# Release Tag Summary for gsDesignNB v0.2.6

## What Was Done

This PR sets up the infrastructure for creating release tags and automated GitHub releases for the gsDesignNB package.

### Files Created/Modified:

1. **`.github/workflows/release.yaml`** - GitHub Actions workflow that:
   - Triggers on push of tags matching `v*.*.*` pattern
   - Extracts version number and release notes from NEWS.md
   - Creates a GitHub release automatically

2. **`create-release-tag.sh`** - Helper script that:
   - Reads current version from DESCRIPTION file
   - Creates properly formatted annotated git tag
   - Optionally pushes tag to GitHub
   - Interactive with user prompts

3. **`RELEASE.md`** - Comprehensive documentation covering:
   - Quick start for v0.2.6 release
   - General release process for future versions
   - CRAN submission workflow
   - Troubleshooting guide

4. **`.Rbuildignore`** - Updated to exclude:
   - RELEASE.md
   - create-release-tag.sh

5. **Local git tag** - Created `v0.2.6` tag on commit 44dde3e

## Current Status

✅ All infrastructure is in place
✅ Local tag `v0.2.6` is created on the correct commit (44dde3e: "Prepare CRAN submission 0.2.6")
✅ Documentation is complete
✅ Helper script is ready to use

⚠️ Tag has NOT been pushed to GitHub yet (cannot be done from this environment)

## Next Steps - Action Required

To complete the release, the repository owner needs to push the tag to GitHub. There are two options:

### Option 1: Automatic (Recommended)
When this PR is merged to master, run:
```bash
./create-release-tag.sh
```

The script will:
1. Detect the v0.2.6 tag already exists locally
2. Ask if you want to recreate it (say No if it's already correct)
3. Push the tag to GitHub
4. The GitHub Actions workflow will automatically create the release

### Option 2: Manual
From the master branch, run:
```bash
git push origin v0.2.6
```

This will push the existing tag and trigger the release workflow.

## Verification

After pushing the tag, you can:
1. Monitor the workflow: https://github.com/keaven/gsDesignNB/actions
2. See the release: https://github.com/keaven/gsDesignNB/releases
3. The release will include notes extracted from NEWS.md

## Important Notes

- The tag points to commit 44dde3e ("Prepare CRAN submission 0.2.6")
- This commit is already prepared for CRAN submission with version 0.2.6
- The NEWS.md already contains appropriate release notes
- No code changes are needed for the release itself

## For Future Releases

For subsequent releases, simply:
1. Update version in DESCRIPTION
2. Add release notes to NEWS.md
3. Run `./create-release-tag.sh`
4. The script handles everything else automatically

---

**Ready to release! Just needs the tag to be pushed to GitHub.**
