# Release Process for gsDesignNB

## Creating a New Release

### Prerequisites
1. Ensure all changes are merged to the master branch
2. Update the `Version` field in `DESCRIPTION` file
3. Update `NEWS.md` with release notes for the new version
4. Run `devtools::check()` to ensure the package passes all checks
5. Commit all changes

### Manual Release Process

1. **Create and push the version tag:**
   ```bash
   # Create an annotated tag with the version number
   git tag -a v0.2.6 -m "Release version 0.2.6 - First CRAN release"
   
   # Push the tag to GitHub
   git push origin v0.2.6
   ```

2. **Automated release creation:**
   - Once the tag is pushed, the `.github/workflows/release.yaml` workflow will automatically:
     - Extract the version number from the tag
     - Extract release notes from `NEWS.md`
     - Create a GitHub release with the tag

### Version Numbering

This package follows [Semantic Versioning](https://semver.org/):
- **Major** (x.0.0): Incompatible API changes
- **Minor** (0.x.0): New functionality in a backward-compatible manner
- **Patch** (0.0.x): Backward-compatible bug fixes

### Current Release: v0.2.6

To create the release for the current version (0.2.6):

```bash
# Tag the commit that prepared CRAN submission
git tag -a v0.2.6 44dde3e -m "Release version 0.2.6 - First CRAN release"

# Push the tag
git push origin v0.2.6
```

This will trigger the release workflow and create a GitHub release with notes extracted from NEWS.md.

### CRAN Submission

After creating the GitHub release:

1. Build the package tarball: `devtools::build()`
2. Run final checks: `devtools::check()`
3. Submit to CRAN via [CRAN submission form](https://cran.r-project.org/submit.html)
4. Respond to any CRAN feedback and update as needed

### Troubleshooting

**If the tag already exists:**
```bash
# Delete the local tag
git tag -d v0.2.6

# Delete the remote tag (if already pushed)
git push origin :refs/tags/v0.2.6

# Recreate the tag
git tag -a v0.2.6 <commit-hash> -m "Release message"
git push origin v0.2.6
```

**If you need to create a tag for a specific commit:**
```bash
git tag -a v0.2.6 <commit-hash> -m "Release message"
git push origin v0.2.6
```
