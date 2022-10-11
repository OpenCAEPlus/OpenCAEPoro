# Automation Scripts

## Directory and Files

- build/:   scripts for building the OpenCAEPoro on MacOS, Windows, and Linux
- publish/: scripts for publishing OpenCAEPoro releases
- doc/:     scripts for building documentations
- other/:   other scripts

## Important

Since our documentation websites requires doxygen process first, then call doxybook2 to convert the xml file to markdown, thus it is easier to manually build the website locally then upload to static site host.

To publish a documentation, you need to follow the steps:

1. `pnpm doc:build:vitepress`
2. commit the changes
3. `pnpm doc:publish:render`