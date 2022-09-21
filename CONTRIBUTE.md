# OpenCAEPoro Contributing Guide

This guide is for developers who want to commit to or create pull requst for this repository.
Please see the README file for user guide.

We are using a trunk based development strategy, which means:
1. Commit to main branch directly when possible, such as when the issue is small and just changing a few lines are enough. Include "closes IssueId" to the commit message. For external developers (cannot commit to main branch), please use the normal github flow, fork and create pull request.
2. When more changes are needed for a functionality, create new branch, but only create short-living branches. Branches should preferrably live less than one day, maybe at most two days. If it is impossible to finish a branch within one day, then you should split your issue into smaller goals, and have the original issue as a project in linear.

Check the [linear doc](https://linear.app/docs/github?tabs=b5eb539099f9#basics) for more on the PR and commit integration with github.

## Commit guide (on main branch)

1. Keep commit small, for example commit every new function added instead of a whole new functionality 
2. As few files as possible to reduce your chance of getting conflicts with others
3. Push multiple times a day to keep the main branch always update and resolve conflicts when it's still managable.
4. Try to keep the HEAD on main branch always buildable
5. Copy issue id from linear and include in your commit message

## Pull Request Guide

1. only have short-living branches
2. user linear generate branch name

## Development Environment Setup

### Code Editor
Recommend using VScode for cross platform compatibility.

### Compiler
Recommend Intel oneapi for easier cross platform compatibility.

GCC setup guide ....

## Documentation Environment Setup

Install texlive ...

Install doxygen...

Install doxygen2...

Install nodejs,npm,pnpm...