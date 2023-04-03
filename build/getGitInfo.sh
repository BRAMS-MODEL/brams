#!/bin/bash
git log -1 > saida.git
if (( $? != 0 ))
then
    echo 'Comando Git nao encontrado'
    echo "!Include para controle de versão" > ../src/utils/include/modGitInfo.inc 
    echo " character(len=*), parameter :: lastCommit  = 'No versioned - did You used  the gitHub  version?'" >> ../src/utils/include/modGitInfo.inc
    echo " character(len=*), parameter :: lastMerge   = 'No versioned - https://github.com/luflarois/brams'" >> ../src/utils/include/modGitInfo.inc
    echo " character(len=*), parameter :: lastAuthor  = 'No versioned - !!!!  Please,  use git clone  !!!!'" >> ../src/utils/include/modGitInfo.inc
    echo " character(len=*), parameter :: lastGitDate = 'No versioned -               ####'" >> ../src/utils/include/modGitInfo.inc   
else
    echo "!Include para controle de versão" > ../src/utils/include/modGitInfo.inc 
    echo " character(len=*), parameter :: lastCommit  = '" $(head -1 saida.git) "'" >> ../src/utils/include/modGitInfo.inc
    echo " character(len=*), parameter :: lastMerge   = '" $(head -2 saida.git | tail -1) "'" >> ../src/utils/include/modGitInfo.inc
    echo " character(len=*), parameter :: lastAuthor  = '" $(head -3 saida.git | tail -1) "'" >> ../src/utils/include/modGitInfo.inc
    echo " character(len=*), parameter :: lastGitDate = '" $(head -4 saida.git | tail -1 ) "'" >> ../src/utils/include/modGitInfo.inc
    rm saida.git
fi