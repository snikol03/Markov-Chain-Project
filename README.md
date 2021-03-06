# Markov Chain

In this project our aim was to use Markov Chain model to classify varius sequences. This project was assigned to us in the course of BIO 331 - Computational and Systems Biology.

_Subject 1:_ Identification of generalized classes of protein molecules.

_Subject 2:_ Identification of proteins with local "abnormal" amino acid composition

_Subject 3:_ Classification of proteins in families

This project was created to check the results of the Markov chain project and compare them with other existing program results.

## Getting Started

You can download this project from the link below:

https://github.com/snikol03/Markov-Chain-Project.git

## Prerequisites

To run this project you need to install perl. You can download and install perl from the following link.

https://www.perl.org/get.html

After the download and installation you can run the program with the terminal or the command line.

#### NOTE:
The .pl file you would like to execute should be in the same folder with the project.pm file.

## Running the tests

To execute this project run the following command in the terminal depending on the perpose you want.

_Subject 1:_ Identification of generalized classes of protein molecules run

`perl thema1.pl`

When this command is executed you will be asked to insert the path of the positive and negative file. Then a menu will appear asking you to select dictionary. At this point there are three options: 1.DNA, 2.RNA and 3.Proteins. The input must be the number of the choice. After that the user have to insert the Markov chain window size. Then an other menu will appear and the user have to choose with a number the tests to run. Depending from the answer the program will need some inputs, except self consistency. For cross validation test the user have to write the number of groups to be devided to. For external data the user have to write the positive and negative external data.

_Subject 2:_ Identification of proteins with local "abnormal" amino acid
composition run

`perl thema2.pl`

When this command is executed you will be asked to insert the path of a .fasta and a .fasta.tab file. Then a menu will appear asking you to select dictionary. At this point there are three options: 1.DNA, 2.RNA and 3.Proteins. The input must be the number of the choice. After that the user have to insert the Markov chain window size. Then an other menu will appear and the user have to choose with a number the tests to run. Depending from the answer the program will need some inputs, except self consistency. For cross validation test the user have to write the number of groups to be devided to. For external data the user have to write the .fasta and .fasta.tab external files.

_Subject 3:_ Classification of proteins in families run

`perl thema4.pl`

When this command is executed you will be asked to insert the path for each family files. When all family files are inserted the user have to write "end" to stop reading. Then a menu will appear asking you to select dictionary. At this point there are three options: 1.DNA, 2.RNA and 3.Proteins. The input must be the number of the choice. After that the user have to insert the Markov chain window size. Then an other menu will appear and the user have to choose with a number the tests to run. Depending from the answer the program will need some inputs, except self consistency. For cross validation test the user have to write the number of groups to be devided to. For external data the user have to write an external file for each family.

For each subject you can select a combination of the following three tests: 
 - Self Consistency
 - Cross Validation
 - External Validation


## Authors

- Alexandra Crysanthou, Mechanical Engineering student, University of Cyprus
- Andria Nicolaou, Electrical Engineering student, University of Cyprus
- Kyriakos Papa, Electrical Engineering student, University of Cyprus
- Sotiroula Afxenti, Computer Science student, University of Cyprus
- Stelios Nikolaou, Computer Science student, University of Cyprus

## License

This project is licensed under the MIT License - see the LICENSE.md file for details
