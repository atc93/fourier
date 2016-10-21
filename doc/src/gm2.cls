%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2009/09/22 v0.2 Internal G-2 Note Class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2014 by Sabato Leo
%
% This file may be distributed and/or modified
%
% 1. under the LaTeX Project Public License and/or
% 2. under the GNU Public License.
%
% See the file doc/licenses/LICENSE for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage example:
%
% \documentclass[letter]{gm2-}
%
% \setlength{\titleblockheight}{6cm}%adjust this lenght if it is needed
%
% \begin{document}
% \title{Template for a G-2 Note}
% \noteid{G-2/Technical/PUBLIC/1}
%
% \version{1.1} %or \draft{0.1}
%
% \author[1,2]{J. Lastname}
% \author[1]{\mail{sleo@illinois.edu}{S. Leo}}
% \author[2]{K. Pitts\thanks{kpitts@illinois.edu}}
% \affil[1]{University Illinois at Champaign-Urbana, IL, USA}
% \affil[2]{University Illinois at Champaign-Urbana, IL, USA}
%
% \maketitle
%
% ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{gm2}[2014/08/10 v0.2 Internal G-2 Note Class]
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[11pt]{article}

\lineskip 1.3pt
\normallineskip 1.3pt
\parskip 2pt
\def\baselinestretch{1.2}\large\normalsize

% REQUIRED PACKAGES
%\RequirePackage[english]{babel}
\RequirePackage{eso-pic,ifpdf,amsmath}
\RequirePackageWithOptions{color,graphicx,forloop}
\RequirePackageWithOptions{hyperref,url,boxedminipage}

% OPTIONS ONLY FOR PDFLaTeX
\ifpdf
 \definecolor{rltred}{rgb}{0.75,0,0}
 \definecolor{rltgreen}{rgb}{0,0.5,0}
 \definecolor{rltblue}{rgb}{0,0,0.75}
 \hypersetup
 {
  bookmarksopen=true,
  colorlinks=true,
  urlcolor=rltblue,
  filecolor=rltgreen,
  linkcolor=rltred
 }
 \pdfimageresolution=300
 \pdfcompresslevel=0
 \newcommand{\mail}[2]{\href{mailto:#1}{#2}}
% OPTIONS ONLY FOR LaTeX (PS OUTPUT)
\else
 \hypersetup
 {
  colorlinks=false%
 }
% \renewcommand{\href}[2]{#2\footnote{#1}}
 \renewcommand{\href}[2]{#2}
 \newcommand{\mail}[2]{\href{#1}{#2}}
\fi

%\pagestyle{headings}

\setlength{\topmargin}{-0.25cm}
\setlength{\textwidth}{15.5cm}
\setlength{\textheight}{22cm}
\setlength{\oddsidemargin}{0.5cm}
\setlength{\evensidemargin}{0.5cm}

\setcounter{topnumber}{9}
\setcounter{bottomnumber}{9}
\setcounter{totalnumber}{9}

\renewcommand{\topfraction}{0.99}
\renewcommand{\bottomfraction}{0.99}
\renewcommand{\textfraction}{0.0}
\renewcommand{\floatpagefraction}{0.99}

\newcommand{\noteid}[1]{\def\@noteid{#1}}

\newcommand{\version}[1]{\def\@version{Version #1}}
\newcommand{\draft}[1]
{
 \def\@version{Draft #1}
 \AddToShipoutPicture{\resizebox{0.9\paperwidth}{0.9\paperheight}%
 {\rotatebox{60}{\color[gray]{0.95}\hspace*{5mm}\textsc{Draft}}}}
}

%    ****************************************
%    *           authblk.sty                *
%    ****************************************
%
%%
%% This is file `authblk.sty',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% authblk.dtx  (with options: `package')
%% =============================================
%% IMPORTANT NOTICE:
%% 
%% This program can be redistributed and/or modified under the terms
%% of the LaTeX Project Public License Distributed from CTAN
%% archives in directory macros/latex/base/lppl.txt; either
%% version 1 of the License, or any later version.
%% 
%% This is a generated file.
%% It may not be distributed without the original source file authblk.dtx.
%% 
%% Full documentation can be obtained by LaTeXing that original file.
%% Only a few abbreviated comments remain here to describe the usage.
%% =============================================
%% Copyright 1997-1999 Patrick W Daly
%% Max-Planck-Institut f\"ur Aeronomie
%% Max-Planck-Str. 2
%% D-37191 Katlenburg-Lindau
%% Germany
%% E-mail: daly@linmpi.mpg.de
%\NeedsTeXFormat{LaTeX2e}[1995/06/01]
%\ProvidesPackage{authblk}
%          [2001/02/27 1.3 (PWD)]
 %-------------------------------------------------------------------
 % NOTICE:
 % This file may be used for non-profit purposes.
 % It may not be distributed in exchange for money,
 %   other than distribution costs.
 %
 % The author provides it `as is' and does not guarantee it in any way.
 %
%% authblk coding copyright (C) 1997 Patrick W. Daly
 % Max-Planck-Institut f\"ur Aeronomie
 % Max-Planck-Str. 2
 % D-37191 Katlenburg-Lindau
 % Germany
 %
 % E-mail: daly@linmpi.mpg.de
 %-----------------------------------------------------------
 % Author/Affiliation Blocks
 % A LaTeX2e package to redefine the \author command to work as normal or
 % to allow a footnote style of author/affiliation input.
 %
 % Three methods of input:
 %  1. Standard LaTeX
 %     \author{Name1//Affil1 /and Name2 and Name3//Affil2 . . .}
 % (use \authorcr, not \\, to force new line within the author text)
 %  2. Automatic mode (selects modes 1 or 3 depending on number of affiliations)
 %     \author{Name1}
 %     \affil{Affil1}
 %     \author{Name2}
 %     \author{Name3}
 %     \affil{Affil2}
 %     . . .
 %  3. Footnotes explicitly
 %     \author[*]{Name1}
 %     \author[**]{Name2}
 %     \author[**]{Name3}
 %     \affil[*]{Affil1}
 %     \affil[**]{Affil2}
 %     . . .
 % See documentation in authblk.dtx for further configuration possibilities.
 %-----------------------------------------------------------
\newcommand\Authfont{\normalfont}
\newcommand\Affilfont{\normalfont}
\newcommand\Authsep{, }
\newcommand\Authands{, and }
\newcommand\Authand{ and }
\newlength{\affilsep}\setlength{\affilsep}{1em}
\newlength{\@affilsep}
\newcounter{Maxaffil}
\setcounter{Maxaffil}{2}
%\DeclareOption{blocks}{\setcounter{Maxaffil}{\maxdimen}}
%\DeclareOption{noblocks}{\setcounter{Maxaffil}{0}}
%\DeclareOption{max2}{\setcounter{Maxaffil}{2}}
%\DeclareOption{max3}{\setcounter{Maxaffil}{3}}
%\DeclareOption{max4}{\setcounter{Maxaffil}{4}}
%\DeclareOption{max5}{\setcounter{Maxaffil}{5}}
%\DeclareOption{max6}{\setcounter{Maxaffil}{6}}
%\DeclareOption{auth-sc}{\renewcommand\Authfont{\scshape}}
%\DeclareOption{auth-sc-lg}{\renewcommand\Authfont{\large\scshape}}
%\DeclareOption{auth-lg}{\renewcommand\Authfont{\large}}
%\DeclareOption{affil-sl}{\renewcommand\Affilfont{\slshape}}
%\DeclareOption{affil-it}{\renewcommand\Affilfont{\itshape}}
%\DeclareOption{german}{\renewcommand\Authand{, und }%
%                       \renewcommand\Authands{, und }}
\ProcessOptions
\newcounter{authors}
\newcounter{affil}
\newif\ifnewaffil \newaffiltrue
\newcommand\AB@authnote[1]{\textsuperscript{\normalfont#1}}
\newcommand\AB@affilnote[1]{\textsuperscript{\normalfont#1}}
\providecommand\textsuperscript[1]{$^{#1}$}
\newcommand\AB@blk@and{\protect\Authfont\protect\AB@setsep}
\newcommand\AB@pand{\protect\and \protect\Authfont \protect\AB@setsep}
\@namedef{@sep1}{}
\@namedef{@sep2}{\Authand}
\newcommand\AB@affilsep{\protect\Affilfont}
\newcommand\AB@affilsepx{\protect\\\protect\Affilfont}
\newcommand\AB@setsep{\setlength{\@affilsep}{\affilsep}}
\newcommand\AB@resetsep{\setlength{\@affilsep}{\z@}}
\newcommand\AB@authlist{}
\newcommand\AB@affillist{}
\newcommand\AB@authors{}
\newcommand\AB@empty{}
\xdef\AB@author{\noexpand\AB@blk@and\@author}
\renewcommand\author[2][]%
      {\ifnewaffil\addtocounter{affil}{1}%
       \edef\AB@thenote{\arabic{affil}}\fi
      \if\relax#1\relax\def\AB@note{\AB@thenote}\else\def\AB@note{#1}%
        \setcounter{Maxaffil}{0}\fi
      \ifnum\value{authors}>1\relax
      \@namedef{@sep\number\c@authors}{\Authsep}\fi
      \addtocounter{authors}{1}%
      \begingroup
          \let\protect\@unexpandable@protect \let\and\AB@pand
          \def\thanks{\protect\thanks}\def\footnote{\protect\footnote}%
         \@temptokena=\expandafter{\AB@authors}%
         {\def\\{\protect\\[\@affilsep]\protect\Affilfont
              \protect\AB@resetsep}%
              \xdef\AB@author{\AB@blk@and#2}%
       \ifnewaffil\gdef\AB@las{}\gdef\AB@lasx{\protect\Authand}\gdef\AB@as{}%
           \xdef\AB@authors{\the\@temptokena\AB@blk@and}%
       \else
          \xdef\AB@authors{\the\@temptokena\AB@as\AB@au@str}%
          \global\let\AB@las\AB@lasx\gdef\AB@lasx{\protect\Authands}%
          \gdef\AB@as{\Authsep}%
       \fi
       \gdef\AB@au@str{#2}}%
         \@temptokena=\expandafter{\AB@authlist}%
         \let\\=\authorcr
         \xdef\AB@authlist{\the\@temptokena
           \protect\@nameuse{@sep\number\c@authors}%
           \protect\Authfont#2\AB@authnote{\AB@note}}%
      \endgroup
      \ifnum\value{authors}>2\relax
      \@namedef{@sep\number\c@authors}{\Authands}\fi
      \newaffilfalse
}
\newcommand\authorcr{\protect\\ \protect\Authfont \protect\AB@setsep}%
\newcommand\affil[2][]%
   {\newaffiltrue\let\AB@blk@and\AB@pand
      \if\relax#1\relax\def\AB@note{\AB@thenote}\else\def\AB@note{#1}%
        \setcounter{Maxaffil}{0}\fi
      \begingroup
        \let\protect\@unexpandable@protect
        \def\thanks{\protect\thanks}\def\footnote{\protect\footnote}%
        \@temptokena=\expandafter{\AB@authors}%
        {\def\\{\protect\\\protect\Affilfont}\xdef\AB@temp{#2}}%
         \xdef\AB@authors{\the\@temptokena\AB@las\AB@au@str
         \protect\\[\affilsep]
         \protect\Affilfont\AB@temp}%
         \gdef\AB@las{}\gdef\AB@au@str{}%
        {\def\\{, \ignorespaces}\xdef\AB@temp{#2}}%
        \@temptokena=\expandafter{\AB@affillist}%
        \xdef\AB@affillist{\the\@temptokena \AB@affilsep
          \AB@affilnote{\AB@note}\protect\Affilfont\AB@temp}%
      \endgroup
       \let\AB@affilsep\AB@affilsepx
}
\def\@author{}
\renewcommand\@author{\ifx\AB@affillist\AB@empty\AB@author\else
      \ifnum\value{affil}>\value{Maxaffil}\def\rlap##1{##1}%
    \AB@authlist\\[\affilsep]\AB@affillist
    \else  \AB@authors\fi\fi}
\let\AB@maketitle=\maketitle
\def\maketitle
  {{\renewenvironment{tabular}[2][]{\begin{center}}
                                   {\end{center}}
  \AB@maketitle}}
%%
%% <<<<< End of generated file <<<<<<
%%
%% End of file `authblk.sty'.
\renewcommand\Affilfont{\fontsize{10}{12}\selectfont}
\renewcommand\AB@affilsepx{; \protect\ \protect\Affilfont}
\setcounter{Maxaffil}{1}
\setlength{\affilsep}{0.3em}

%% A slightly modified definition of \maketitle
%%
\newlength{\titleblockheight}       % so user can change it if need be
\setlength{\titleblockheight}{5.5cm}

\def\@maketitle{\vbox to \titleblockheight {\hsize\textwidth
  \linewidth\hsize \vskip 3pt
  \begin{flushright}
   {\large\tt\@noteid}\\
   \@version\\
   \today
  \end{flushright}
  \vskip 1em
  \centering
  {\Large\bf \@title \par}
  \vskip 2em                % Vertical space after title.
  {\large\begin{tabular}[t]{@{}c@{}}\@author \end{tabular}\par}
  \vfil}}