#!/bin/bash

# SHELL script for building the documentation and push to github pages
# taken from here: https://www.alkaline-ml.com/2018-12-23-automate-gh-builds/

set -e
shopt -s extglob

cd doc
make clean html
cd ..

mv doc/_build/html ./
rm -r !(".git"|"html"|".."|".")

git checkout gh-pages
git checkout --orphan gh-pages-tmp
git config --global user.email "$GH_EMAIL" > /dev/null 2>&1
git config --global user.name "$GH_NAME" > /dev/null 2>&1
touch .nojekyll

if [[ "$CIRCLE_BRANCH" =~ ^master$|^[0-9]+\.[0-9]+\.X$ ]]; then
    cp -r html/* ./
else
    mkdir -p "$CIRCLE_BRANCH"
    cp -r html/* ./"$CIRCLE_BRANCH"
fi

rm -r html/    
git add --all
git commit --allow-empty  -m "[ci skip] Publishing updated documentation..."
git remote rm origin
git remote add origin https://"$GH_NAME":"$GH_TOKEN"@github.com/yqshao/mdppp.git
git push --force origin gh-pages-tmp:gh-pages
