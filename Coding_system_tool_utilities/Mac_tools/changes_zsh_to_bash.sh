https://crunchify.com/how-to-change-mac-os-x-terminal-color/

Open terminal window
Click on Terminal Menu
Click on Preference
Modify profile as per below image
Select profile for



On startup, 

1) select "text" change the color as you like  

2) selet the "shell" to run 
source .bashrc


# mac change zsh to bash
# https://medium.com/macoclock/how-to-change-the-colour-of-your-bash-prompt-on-mac-b06032543353
chsh -s /bin/bash
export PS1="\h:\W \[\033[1;31m\]\u\\[\033[0m\]$ "


# https://apple.stackexchange.com/questions/33677/how-can-i-configure-mac-terminal-to-have-color-ls-output

export CLICOLOR=1  

alias l="ls -CF"
alias ll="ls -l"
alias la="ls -a"
alias lla="ls -la"
alias dir='dir --color=auto'
alias vdir='vdir --color=auto'
alias grep='grep --color=auto'
alias fgrep='fgrep --color=auto'
alias egrep='egrep --color=auto'
# some more ls aliases
alias c='clear'
alias h='history'


GREEN='\e[1;32m'
YELLOW='\e[33m'
BLUE='\e[1;34m'
PURPLE='\e[1;35m'
CYAN='\e[1;36m'
RESET="\[$(tput sgr0)\]"

color_prompt=yes

new_line() {
    printf "\n$ "
}
# choose any one the the following 3 
PS1="\h:\W \[\033[1;31m\]\u\\[\033[0m\]$ "
PS1="\w : \[\033[1;31m\]\u:\h \n$\[\033[0m\] " # : can be replaced with @ 
PS1="\[\033[1;38m\]\w : \[\033[0m\] \n\[\033[1;31m\]\u@\h $\[\033[0m\] "
PS1="\[\035[1;35m\]\w : \[\033[0m\] \n\[\033[1;31m\]\u:\h $\[\033[0m\] " # magnet  or prupple
PS1="\[\032[1;32m\]\w : \[\033[0m\] \n\[\033[1;31m\]\u:\h $\[\033[0m\] " # Green

PS1="\[\032[1;92m\]\w : \[\033[0m\] \n\[\033[1;31m\]\u:\h $\[\033[0m\] " # light green
PS1="\[\032[1;96m\]\w : \[\033[0m\] \n\[\033[1;31m\]\u:\h $\[\033[0m\] " # light cyan

PS1="\[\032[1;42m\]\w : \[\033[0m\] \n\[\033[1;31m\]\u:\h $\[\033[0m\] " # with backgroud Green
PS1="\[\032[1;47m\]\w : \[\033[0m\] \n\[\033[1;31m\]\u:\h $\[\033[0m\] " # with backgroud Grey

# green dir and red next line
PS1="\[\033[1;92m\]\w : \[\033[0m\] \n\[\033[1;31m\]\u@\h $\[\033[0m\] "
### good 1) purple pwd, 2) green User, 3) red Hoster
PS1="\[\033[1;35m\]\w : \[\033[0m\] \n\[\033[1;32m\]\u@\[\033[0m\] \[\033[1;31m\]\h $\[\033[0m\] "
### good 1) purple pwd, 2) blue User, 3) red Hoster
PS1="\[\033[1;35m\]\w : \[\033[0m\] \n\[\033[1;34m\]\u@\[\033[0m\] \[\033[1;31m\]\h $\[\033[0m\] "

#https://misc.flogisoft.com/bash/tip_colors_and_formatting

#https://www.makeuseof.com/customize-zsh-prompt-macos-terminal/
#It looks all setting was in the 2 file : 
/etc/zshrc
/etc/zshrc_Apple_Terminal


+++=No previlige to change it

