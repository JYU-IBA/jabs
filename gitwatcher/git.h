#ifndef JABS_GIT_H
#define JABS_GIT_H

int git_populated();
int git_dirty();
const char *git_commit_sha1();
const char *git_commit_date();
const char *git_describe();
const char *git_branch();
#endif //JABS_GIT_H
