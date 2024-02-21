#ifndef JABS_GIT_H
#define JABS_GIT_H

int git_populated(void);
int git_dirty(void);
const char *git_commit_sha1(void);
const char *git_commit_date(void);
const char *git_describe(void);
const char *git_branch(void);
#endif //JABS_GIT_H
