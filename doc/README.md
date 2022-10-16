# Documentation

## Directory and Files

- doxygen, folder for .dox file complimenting the documentation in source code
- doxygen-site, folder for holding the doxygen generated website, [https://opencaeplus.github.io/OpenCAEPoro/](https://opencaeplus.github.io/OpenCAEPoro/)
- website, contain the markdown files and vitepress configuration
- Doxyfile, doxygen configuration file
- doxybook-config.json, doxybook2 configuration file, for converting doxygen xml to markdown
- header.tmpl, header template for doxybook2

## Write

### Doxygen

### Markdown

## Build

To build doc, you need to first setup the development environment

### Environment Setup

1. doxygen: for parsing the source code
2. latex: for generating pdf from doxygen generated tex file
3. doxybook2: for creating markdown file from doxygen xml
4. pnpm: for creating the website
