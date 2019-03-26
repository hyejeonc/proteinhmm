## Hiddem Markov Model : Prediction of protein secondary structure 
##### This is a repository for project 'Prediction of protein secondary structure from sequence by Hiddem Markov Model(HMM) method' for a course 'MOL3022- Bioinformatics - Method Oriented Project' at NTNU, in spring 2019. 
There are many folders and files, but only **'final'** folder is related to the complete program.


##### How to install 
- hiddenmarkovmodel.py  -- main file 
- statistics.py-- module file 
- sequence.py -- module file
- mathematics.py -- module file
  
###### Windows
1. Download all files in /final/windows/ folder. 
2. Excecute hiddenmarkovmodel.exe file. 
3. Type file directory for training and testing.
4. Output files (decoded sequence) are saved at output folder.
(Caution! 'output' folder MUST exist for saving.)

###### Linux
1. Download four files above at ./final/linux/ 
2. Input files MUST be in same folder with hiddenmarkovmodel.py 
3. Type terminal command below at the directory that you have four .py files mentioned above.
```
chmod u+x hiddenmarkovmodel.py
```
3. Type terminal command below for executing.
```
./hiddenmarkovmodel.py
```
4. Output files are saved in same directory with hiddenmarkovmodel.py


##### Main files(./final/)
- [x] hiddenmarkovmodel.py  -- main file -- last update 25 Mar, 2019
- [ ] statistics.py-- module file -- last update 22 Mar, 2019 
- [x] sequence.py -- module file -- last update 25 Mar, 2019
- [x] mathematics.py -- module file -- last update 25 Mar, 2019

##### Datasets(./final/dataset/)
- protein-secondary-structure.train -- for training, 111 sets
- protein-secondary-structure.test -- for testing, 17 sets

##### Outputs(./final/output/)
- raw_data_simpleHMM.csv
- raw_data_EMHMMconv.csv
- raw_data_EMHMMiter.csv

##### Reference documents and links
1. Jones, D.T.: Protein secondary structure prediction based on position-specic scoring matrices. Journal of
molecular biology 292(2), 195{202 (1999)
2. Example Diagram of HMM. http://yanfenglu.net/researchVAS_p1.htm
3. Abela, J., Michael, J.: Topics in evolving transformation systems [microform]. (2019)
4. Sheh, A., Ellis, D.P.: Chord segmentation and recognition using em-trained hidden markov models (2003)
5. Bayes Theorem. https://www.britannica.com/topic/Bayess-theorem
6. Hidden Markov Models.
https://archive.ics.uci.edu/ml/datasets/Molecular+Biology+(Protein+Secondary+Structure)
7. Qian, N., Sejnowski, T.J.: Predicting the secondary structure of globular proteins using neural network models.
Journal of molecular biology 202(4), 865{884 (1988)
8. Main Program Repository of Author. https://github.com/hyejeonc
9. Matlib Library. https://matplotlib.org/
10. Asai, K., Hayamizu, S., Handa, K.: Prediction of protein secondary structure by the hidden markov model.
Bioinformatics 9(2), 141{146 (1993)
11. Boodidhi, S.: Using smoothing techniques to improve the performance of hidden markov's model (2011)
12. Shifrin, J., Pardo, B., Meek, C., Birmingham, W.: Hmm-based musical query retrieval. In: Proceedings of the
2nd ACM/IEEE-CS Joint Conference on Digital Libraries, pp. 295{300 (2002). ACM
13. Landschulz, W.H., Johnson, P.F., McKnight, S.L.: The leucine zipper: a hypothetical structure common to a
new class of dna binding proteins. Science 240(4860), 1759{1764 (1988)
14. Soding, J.: Protein homology detection by hmm{hmm comparison. Bioinformatics 21(7), 951{960 (2004)
15. Won, K.-J., Hamelryck, T., Prugel-Bennett, A., Krogh, A.: An evolutionary method for learning hmm structure:
prediction of protein secondary structure. BMC bioinformatics 8(1), 357 (2007)
16. Lee, L., Leopold, J.L., Frank, R.L.: Protein secondary structure prediction using blast and relaxed threshold rule
induction from coverings. In: 2011 IEEE Symposium on Computational Intelligence in Bioinformatics and
Computational Biology (CIBCB), pp. 1{8 (2011). IEEE
17. Protein HMM raw data files by programming tool https://github.com/dmvvliet/protein-HMMs
18. hmmexample.py -- main file for predicting weather with comsumption of icecream https://github.com/jason2506/PythonHMM
19. weatherexample.py -- module file, same as 18.
20. COMP3212CompBiologyLabs-master -- predicting genetic code from protein sequence https://github.com/aloisklink/COMP3212CompBiologyLabs
21. HMM and EM algorithm lecture pdf https://people.cs.umass.edu/~mccallum/courses/inlp2004a/lect10-hmm2.pdf
22. HMM with weather example and tutorial https://sambaiga.github.io/2017/05/03/hmm-intro.html
23. Baum-Welch(EM) algorithm http://www.biostat.jhsph.edu/bstcourse/bio638/notes/HMMs_BaumWelch.pdf
24. Smoothing method https://danieltakeshi.github.io/2015-07-25-hidden-markov-models-and-particle-filtering/
