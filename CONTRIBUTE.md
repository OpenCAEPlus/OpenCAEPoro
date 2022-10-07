# OpenCAEPoro Contributing Guide

This guide is for developers who wish to commit to or create pull requests for this repository.

Please refer to the README.md file for a simple user guide.

We are using a `trunk-based development strategy`, which means:

1. Commit to the main branch directly when possible, such as when the issue is small and just changing a few lines are enough. Include "closes IssueId" in the commit message. For external developers (who cannot commit to the main branch directly), please use the normal github workflow--fork and create pull request.

2. When more changes are needed for a functionality, create a new branch, but only create short-living branches. Branches should preferably live less than one or two days. If it is impossible to finish a branch within one day, then you should split your issue into smaller goals, and have the original issue as a project in linear.

Check the [Linear Doc](https://linear.app/docs/github?tabs=b5eb539099f9#basics) for more details on the PR and commit integration with github.

## Commit Guide 

When commit to the main branch, try to
1. Keep commits small, for example commit every new function added instead of a whole new functionality;
2. Change as few files as possible to reduce your chance of getting conflicts with others;
3. Push multiple times a day to keep the main branch always update and resolve conflicts when it's still manageable;
4. Keep the HEAD on the main branch always buildable;
5. Copy IssueId from Linear and include it in your commit message.

## Pull Request Guide

1. Only keep short-living branches;
2. Use [Linear](https://linear.app/) generate branch name.

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

GCC setup guide ....

### CMake

Make sure you adopt the modern cmake approach when writing CMakeLists.txt file. Here is a list of good articles you should read, [modern cmake](https://cliutils.gitlab.io/modern-cmake/), [some modern cmake tips](https://www.incredibuild.com/blog/modern-cmake-tips-and-tricks), [on target_sources](https://crascit.com/2016/01/31/enhanced-source-file-handling-with-target_sources/), [on PRIVATE PUBLIC INTERFACE](https://kubasejdak.com/modern-cmake-is-like-inheritance)

## Documentation Environment Setup

1. Install TeXLive ...

2. Install doxygen ...

3. Install doxygen2 ...

4. Install nodejs,npm,pnpm ...