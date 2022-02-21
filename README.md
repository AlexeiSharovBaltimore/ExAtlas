# ExAtlas
Analysis of gene expression and gene set data
Developed at the National Institute on Aging, NIA/NIH, Baltimore USA.
Programmer Alexei Sharov sharov@comcast.net

This software is provided "AS IS".  Programmer makes no warranties, express or implied, including no representation or warranty with respect to
the performance of the software and derivatives or their safety, effectiveness, or commercial viability. Programmer does not warrant the
merchantability or fitness of the software and derivatives for any particular purpose, or that they may be exploited without infringing
the copyrights, patent rights or property rights of others. Programmer shall not be liable for any claim, demand or action for any loss, harm,
illness or other damage or injury arising from access to or use of the software or associated information, including without limitation any
direct, indirect, incidental, exemplary, special or consequential damages. This software program may not be sold, leased, transferred, exported
or otherwise disclaimed to anyone, in whole or in part, without the prior written consent of programmer.

DIRECTORY STRUSTURE

ExAtlas = exatlas.ini file

 |
 
 ----exatlas = html files
 
 |      |
 
 |      ----bin = exatlas.cgi (CGI program)
 
 |      |
 
 |      ----images = images
 
 |      |
 
 |      ----output = writable directory for output files
 
 |      |
 
 |      ----download
 
 |      
 
 |     
 
 ----exatlasInfo
 
        |
        
        ----bin = other programs (C, perl), C-programs should be compiled here
        
        |
        
        ----info = writable directory for personal config files, login.txt
        
        |
        
        ----data = writable directory for data files
        
        

Compile C-programs in ExAtlas/exatlasInfo/bin:
gcc anova_oneway.c -lm -o anova_oneway
gcc togif.c -lm -o togif
gcc correlation_exatlas.c -lm -o correlation_exatlas
gcc norm_new.c -lm -o norm_new
gcc page_exatlas.c -lm -o page_exatlas
gcc pairwise.c -lm -o pairwise
gcc pca.c -lm -o pca

Configure two accounts to your name: administrator, public
Edit file login.txt and put your name and email address
Log into ExAtlas and change your passwords for both account
Create your own working account in ExAtlas
