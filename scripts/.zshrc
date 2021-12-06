# If you come from bash you might have to change your $PATH.
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/3seq:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/structure:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/tomahawk:$PATH"
export PATH="/Applications/BEAST_2.6.0/bin:$PATH"
export PATH=$HOME/bin:/usr/local/bin:$PATH
export PATH="/Users/ChrisOwen/anaconda/bin:$PATH"
export PATH="/usr/local/opt/gnuplot@4/bin:$PATH"
export PATH="/usr/local/Cellar/tree/1.8.0/bin:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/iqtree-1.6.9-MacOSX/bin:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/IGV_2.4.16:$PATH"
export PATH="/usr/local/Cellar/qt/5.12.0/bin:$PATH"
#export PATH="/usr/local/Cellar/perl/5.28.1/bin:$PATH"
export PATH="/usr/local/Cellar/python/3.6.5/bin:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/quast-5.0.2:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/ratt-code:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/MeShClust-master/bin:$PATH"
export PATH="/usr/local/Cellar/capnp/0.7.0/bin:$PATH"
export PATH="/usr/local/Cellar/autoconf/2.69/bin:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/Mash-master:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/EMBOSS-6.6.0/emboss:$PATH"
export PATH="/Applications/CMake.app/Contents/bin:$PATH"
export PATH="/usr/local/opt/llvm/bin:$PATH"
export PATH="/usr/local/opt/icu4c/bin:$PATH"
export PATH="/usr/local/opt/icu4c/sbin:$PATH"
export PATH="/usr/local/Cellar/hdf5/1.10.5_1/bin:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/bioawk-master:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/bindash-master/release:$PATH"
export PATH="/usr/local/ncbi/blast/bin:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/FastTree:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/trimal-trimAl/source:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/plink_mac_20200616:$PATH"
export PATH="/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/software/admixture_macosx-1.3.0:$PATH"
source ~/perl5/perlbrew/etc/bashrc
export PERL5LIB="$PERL5LIB:$HOME/lib/perl5"

export TERM=xterm-256color

# Path to your oh-my-zsh installation.
export ZSH="/Users/ChrisOwen/.oh-my-zsh"

# Set name of the theme to load --- if set to "random", it will
# load a random theme each time oh-my-zsh is loaded, in which case,
# to know which specific one was loaded, run: echo $RANDOM_THEME
# See https://github.com/robbyrussell/oh-my-zsh/wiki/Themes
ZSH_THEME="bira"

# Set list of themes to pick from when loading at random
# Setting this variable when ZSH_THEME=random will cause zsh to load
# a theme from this variable instead of looking in ~/.oh-my-zsh/themes/
# If set to an empty array, this variable will have no effect.
# ZSH_THEME_RANDOM_CANDIDATES=( "robbyrussell" "agnoster" )

# Uncomment the following line to use case-sensitive completion.
# CASE_SENSITIVE="true"

# Uncomment the following line to use hyphen-insensitive completion.
# Case-sensitive completion must be off. _ and - will be interchangeable.
# HYPHEN_INSENSITIVE="true"

# Uncomment the following line to disable bi-weekly auto-update checks.
# DISABLE_AUTO_UPDATE="true"

# Uncomment the following line to change how often to auto-update (in days).
# export UPDATE_ZSH_DAYS=13

# Uncomment the following line to disable colors in ls.
# DISABLE_LS_COLORS="true"

# Uncomment the following line to disable auto-setting terminal title.
# DISABLE_AUTO_TITLE="true"

# Uncomment the following line to enable command auto-correction.
# ENABLE_CORRECTION="true"

# Uncomment the following line to display red dots whilst waiting for completion.
COMPLETION_WAITING_DOTS="true"

# Uncomment the following line if you want to disable marking untracked files
# under VCS as dirty. This makes repository status check for large repositories
# much, much faster.
# DISABLE_UNTRACKED_FILES_DIRTY="true"

# Uncomment the following line if you want to change the command execution time
# stamp shown in the history command output.
# You can set one of the optional three formats:
# "mm/dd/yyyy"|"dd.mm.yyyy"|"yyyy-mm-dd"
# or set a custom format using the strftime function format specifications,
# see 'man strftime' for details.
HIST_STAMPS="dd.mm.yyyy"

# Turn off all beeps
unsetopt BEEP
# Turn off autocomplete beeps
# unsetopt LIST_BEEP

# Would you like to use another custom folder than $ZSH/custom?
# ZSH_CUSTOM=/path/to/new-custom-folder

# Which plugins would you like to load?
# Standard plugins can be found in ~/.oh-my-zsh/plugins/*
# Custom plugins may be added to ~/.oh-my-zsh/custom/plugins/
# Example format: plugins=(rails git textmate ruby lighthouse)
# Add wisely, as too many plugins slow down shell startup.
plugins=(
	git
	zsh-syntax-highlighting
	zsh-autosuggestions
)

source $ZSH/oh-my-zsh.sh

# User configuration

# export MANPATH="/usr/local/man:$MANPATH"

# You may need to manually set your language environment
# export LANG=en_US.UTF-8

# Preferred editor for local and remote sessions
# if [[ -n $SSH_CONNECTION ]]; then
#   export EDITOR='vim'
# else
#   export EDITOR='mvim'
# fi

# Compilation flags
# export ARCHFLAGS="-arch x86_64"

# ssh
# export SSH_KEY_PATH="~/.ssh/rsa_id"

# Set personal aliases, overriding those provided by oh-my-zsh libs,
# plugins, and themes. Aliases can be placed here, though oh-my-zsh
# users are encouraged to define aliases within the ZSH_CUSTOM folder.
# For a full list of active aliases, run `alias`.
#
# Example aliases
# alias zshconfig="mate ~/.zshrc"
# alias ohmyzsh="mate ~/.oh-my-zsh"
#alias gcc="/usr/local/Cellar/gcc/11.2.0/bin/gcc-11"
#alias g++="/usr/local/Cellar/gcc/11.2.0/bin/g++-11"
#alias c++="/usr/local/Cellar/gcc/11.2.0/bin/c++-11"
alias tails="ssh chriowen@tails.cs.ucl.ac.uk"
alias comic="ssh chriowen@comic.cs.ucl.ac.uk"
alias bchuckle="ssh chriowen@bchuckle.cs.ucl.ac.uk"
alias pchuckle="ssh chriowen@pchuckle.cs.ucl.ac.uk"
alias j11="export JAVA_HOME=`/usr/libexec/java_home -v 11`; java -version"
alias j8="export JAVA_HOME=`/usr/libexec/java_home -v 1.8`; java -version"
alias open='open -a Finder ./'              # f: Opens current directory in MacOS Finder
#alias cp='cp -iv'                           # Preferred 'cp' implementation
#alias mv='mv -iv'                           # Preferred 'mv' implementation
alias mkdir='mkdir -pv'                     # Preferred 'mkdir' implementation
alias up='cd ../'                           # Go back 1 directory level (for fast typers)
alias upp='cd ../../'                       # Go back 2 directory levels
alias uppp='cd ../../../'                   # Go back 3 directory levels
alias upppp='cd ../../../../'               # Go back 4 directory levels
alias uppppp='cd ../../../../../'           # Go back 5 directory levels
alias upppppp='cd ../../../../../../'       # Go back 6 directory levels
alias find="find . -name "                  # qfind: Quickly search for file
alias fasta2tbl='/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/scripts/fasta2tbl.txt'
alias tbl2fasta='/Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/scripts/tbl2fasta.txt'
alias gbextract='python3 /Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/scripts/gb_extract.py'
alias envall='( setopt posixbuiltin; set ; ) | less'

# directory shortcuts
alias cauli='cd /Users/ChrisOwen/Dropbox/Work/PhD/Bioinfo/data/plant_viruses/Caulimoviridae' 
alias rvh='cd /Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus'

# Custom Commands
# cdl: change directory and list what's inside
cdl () { cd "$1" && ls ; }
# mcd: Makes new Dir and jumps inside
mcd () { mkdir -p "$1" && cd "$1" ; }
# trash: Moves a file to the MacOS trash
trash () { command mv "$@" ~/.Trash ; }
# stripn: makes 2 line fasta format
stripn () { sed -e 's/\(^>.*$\)/#\1#/' "$1" | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' ; }
# seqlen: prints individual sequnece lengths of a multifasta
seqlen () { bioawk -c fastx '{ print $name, length($seq) }' < "$1" }
# splitfasta: split a multifasta into individual seqs with seqname as filename (2 line format)
splitfasta () { cat "$1" | awk 'BEGIN{FS=" ";RS=">"}{ filename=($1".fasta");print ">" $0 > filename;close(filename)}' }
# keep_ids: prune a multifasta based on a list of seq IDs to keep (can be fuzzy match) 2 line format
keep_ids () { fasta2tbl "$1" | grep -f "$2" | tbl2fasta | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' ; }	
# drop_ids: same as above, but invert match, so supply list of seqs to exclude. NOTE: list of patterns supplied MUST NOT contain empty lines - reults in empty file
drop_ids () { fasta2tbl "$1" | grep -vf "$2" | tbl2fasta | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' ; }
eval "$(jenv init -)"

##>>> conda initialize >>>
## !! Contents within this block are managed by 'conda init' !!
#__conda_setup="$('/Users/ChrisOwen/anaconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    eval "$__conda_setup"
#else
#    if [ -f "/Users/ChrisOwen/anaconda3/etc/profile.d/conda.sh" ]; then
#        . "/Users/ChrisOwen/anaconda3/etc/profile.d/conda.sh"
#    else
#        export PATH="/Users/ChrisOwen/anaconda3/bin:$PATH"
#    fi
#fi
#unset __conda_setup
##<<< conda initialize <<<

export PATH="/usr/local/opt/curl/bin:$PATH"
