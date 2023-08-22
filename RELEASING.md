# Releasing level1c4pps

prerequisites: `pip install setuptools twine`


1. checkout main
2. pull from repo
3. run the unittests
4. run `loghub` and update the `CHANGELOG.md` file:

```
loghub foua-pps/level1c4pps --token <personal access token (see https://github.com/settings/tokens)>  -st v<previous version> -plg bug "Bugs fixed" -plg enhancement "Features added" -plg documentation "Documentation changes"
```

Example:
```
loghub foua-pps/level1c4pps -st v0.2.28 -plg bug "Bugs fixed" -plg enhancement "Features added" -plg documentation "Documentation changes"
```

Or only:
```
loghub foua-pps/level1c4pps  -st v0.2.31 
```

The  "--token <personal access token (see https://github.com/settings/tokens)>" somtimes can be excluded.

Paste CHANGELOG.temp to CHANGELOG.md and update RELEASE-VERSION to the version that will be released.

Don't forget to commit!

5. Create a tag with the new version number, starting with a 'v', eg:

```
git tag -a v0.2.29 -m "Version 0.2.29"
```

See [semver.org](http://semver.org/) on how to write a version number.


6. push changes to github `git push name-of-foua-pps-repo main --follow-tags`

7. Verify the Github actions unit tests passed

8. Create a "Release" on GitHub by going to
   https://github.com/foua-pps/level1c4pps/releases and clicking "Draft a new
   release". On the next page enter the newly created tag in the "Tag version"
   field, "Version X.Y.Z" in the "Release title" field, and paste the markdown
   from the changelog (the portion under the version section header) in the
   "Describe this release" box. Finally click "Publish release".

9. Verify the GitHub actions for deployment succeed and the release is on PyPI
