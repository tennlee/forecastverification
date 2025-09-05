# Contributing Guide

> **Note:** This contributor's guide is currently a work in progress.

This repository contains the HTML-based WMO JWGFVR (Joint Working Group on Forecast Verification Research) verification webpage, previously hosted at [https://www.cawcr.gov.au/projects/verification/](https://www.cawcr.gov.au/projects/verification/). This document describes the process for contributors who want to make changes to the website.

## Getting Started
### Fork & Clone
1. Fork the repository on GitHub.
2. Clone your fork:
   ```bash
   git clone https://github.com/<your-github-username>/forecastverification.git
   cd forecastverification
   ```
3. Add upstream
    ```bash
    git remote add upstream https://github.com/JWGFVR/forecastverification.git
    ```
## Preview Locally
1. Open `index.html` or, 
2. Run `python3 -m http.server` and navigate to `http://localhost:8000`

## Contribution Workflow
1. Create a GitHub issue describing the change you want to make:
   - Open an issue at: https://github.com/JWGFVR/forecastverification/issues
   - If an issue already exists for your change, mention it in your PR.
2. Sync your fork with upstream (recommended before creating a branch):
   ```bash
   # fetch and update your local main
   git fetch upstream
   git checkout main
   git pull upstream main
3. **Branch**: Create a branch with a descriptive name
   ```bash
   git checkout -b issue-123-add-crps-score
   ```
4. Develop: Make enhancements to the website.
5. Preview: Test changes locally before committing. 
6. Stage and commit your changes with a clear commit message. E.g.,
    ```bash
    git add index.html
    git commit -m "feat: Added CRPS equation"
    ```
7. Push your branch to your fork
    ```bash
    git push -u origin <branch-name>
    ```
8. On GitHub, open a pull request from your-username:issue-123-add-crps → JWGFVR:main. In the pull request description, reference the issue (e.g., Closes \#123) and explain what you changed and why.
9. Address review feedback:
    - Respond to comments and push updates to the same branch.
    - Keep your branch up to date with upstream main while the PR is open:
    ```bash
    # fetch upstream
    git fetch upstream
    # rebase onto latest main
    git checkout issue-123-add-crps
    git rebase upstream/main
    # resolve conflicts if any, then
    git push --force-with-lease
    ```
10. Merge: After at least one maintainer approves, a maintainer will merge the pull request.

## Content guidelines
TODO. Add some guidelines around what type of content is suitable to add.
