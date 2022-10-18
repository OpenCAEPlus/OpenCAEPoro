# OpenCAEPoro Contributing Guide

This guide is for developers who wish to commit to or create pull requests for this repository.

Please refer to the README.md file for a simple user guide.

## General Guide

### Workflow Explained

1. Create an issue and discuss with others how to fix a bug or add a feature you need;
2. Either fork (external developer) or branch (internal developer) to develope new code and keep your branch life span as short as possible;
3. Create a pull request;
4. An action will automatically create a release draft;
5. Once steps 1-4 are repeated multiple times and a new version is ready, repo admin will publish the drafted release and a changelog will be generated automatically.

### Commit Guide

1. Keep commits small, for example commit every new function added instead of a whole new functionality;
2. Change as few files as possible to reduce your chance of getting conflicts with others;
3. Push multiple times a day to keep the main branch always update and resolve conflicts when it's still manageable;
4. Keep the HEAD on the main branch always build-able;
5. Copy IssueId from Linear (if you are using `Linear`) and include it in your commit message.

### Pull Request Guide

1. Only keep short-living branches/forks; Branches should preferably live less than one or two days. If it is impossible to finish a branch within one day, then you should split your issue into smaller goals.

## For Internal Developers

We are using a mixture of `trunk-based development strategy` and the `github flow`, which means:

1. For anything that you would like to appear on the release log and change log, you should create a short living branch and pull request;
2. For anything that doesn't need to appear in change log, you can commit to the main branch directly, such as when the issue is small, nobody has reported it yet, and just changing a few lines are enough.

<!-- Check the [Linear Doc](https://linear.app/docs/github?tabs=b5eb539099f9#basics) for more details on the PR and commit integration with github. -->

## For External Developers

Use the `github flow`, which means:

1. Fork the repo;
2. Develope the code and pull request

## Development Environment Setup

### Code Style

We list a few rule of thumbs for styling the code here:

1. Use .clang-format to format the code automatically;
2. Naming conventions:
    - Files, classes, methods, functions: Upper Camel;
    - Variables: camel;
    - Constants: all upper case;
    - Make your names meaningful.
3. Keep all functions short and make them do only one thing;
4. Provide Doxygen-style comments.

### Code Editor

Recommend using VScode for cross-platform compatibility.

### Compiler

Recommend Intel oneAPI for easier cross-platform compatibility.

### CMake

Make sure you adopt the modern cmake approach when writing CMakeLists.txt file. Here is a list of good articles you should read, [modern cmake](https://cliutils.gitlab.io/modern-cmake/), [some modern cmake tips](https://www.incredibuild.com/blog/modern-cmake-tips-and-tricks), [on target_sources](https://crascit.com/2016/01/31/enhanced-source-file-handling-with-target_sources/), [on PRIVATE PUBLIC INTERFACE](https://kubasejdak.com/modern-cmake-is-like-inheritance)

## Documentation Environment Setup

1. Install TeXLive ... (TBA)

2. Install doxygen ... (TBA)

3. Install doxygen2 ... (TBA)

4. Install nodejs,npm,pnpm ... (TBA)
