#!/bin/bash
# Script to create and push a release tag for gsDesignNB
# Usage: ./create-release-tag.sh

set -e

# Get current version from DESCRIPTION file
VERSION=$(grep "^Version:" DESCRIPTION | sed 's/Version: //')

echo "Current version in DESCRIPTION: $VERSION"
echo ""

# Check if tag already exists
if git rev-parse "v$VERSION" >/dev/null 2>&1; then
    echo "Tag v$VERSION already exists."
    read -p "Do you want to delete and recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Deleting existing tag v$VERSION..."
        git tag -d "v$VERSION"
        git push origin ":refs/tags/v$VERSION" 2>/dev/null || echo "Tag not found on remote"
    else
        echo "Aborted."
        exit 1
    fi
fi

# Create annotated tag
echo "Creating annotated tag v$VERSION..."
RELEASE_NOTE=$(sed -n "/^# gsDesignNB $VERSION/,/^# gsDesignNB /p" NEWS.md | sed '$d' | tail -n +3 | head -1 | sed 's/^- //')
git tag -a "v$VERSION" -m "Release version $VERSION - $RELEASE_NOTE"

echo ""
echo "Tag v$VERSION created successfully!"
echo ""
echo "To push the tag and trigger the release workflow:"
echo "  git push origin v$VERSION"
echo ""
read -p "Push the tag now? (y/N): " -n 1 -r
echo

if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Pushing tag v$VERSION..."
    git push origin "v$VERSION"
    echo ""
    echo "Tag pushed! The release workflow should now run automatically."
    echo "Check: https://github.com/keaven/gsDesignNB/actions"
else
    echo "Tag created locally but not pushed."
    echo "Push it later with: git push origin v$VERSION"
fi
