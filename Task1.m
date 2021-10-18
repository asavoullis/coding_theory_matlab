clear all; close all; clc
%%TO SHOW ANSWERS REMOVE THE ; FROM THE END OF THE FUNCTIONS

%%TASK1 (15,11) Hamming code

%% Task 1.1: Construction of G and H matrices

% (n,k) -> (15,11)
%k = length of the input sequence.
kin = 11;
%n = length of the output of the encoder.
nout = 15;

%Generates the G and H matrices
[G,H] = hamming_15_11_G_H_matrices();


%%  TASK 1.2 Encoding

%Input = uob number 1637541 / 0805547 / 1145506
UoB_ID = 1637541;

%Encodes the University of Bristol ID
Codeword = UOB_encoding_HC_UoBID(UoB_ID,G)


%% Task 1.3 Decoding

%Constructs the table and the syndrome
[table,syndrome] = table_syndrome_construction(nout,kin,H);

%Decode the 
[Codeword_out,Fixcodeword,CorrectCW] = decoding_codeword(Codeword,UoB_ID,H,table,syndrome);

CorrectCW; %checking

%% Task 1.4 Performance of (15,11) Hamming Code

%Number of Samples 
N = 20000;
%make a boolean to change input you want to use
checkifcorrect = true;
%create a matrix with the specifications:
%only use 1 or 2,rows/samples(N),columns(11),-1 to convert to binary
%this creates an array with various different codewords.
input = randi(2,N,11,'double')-1;

%generate n samples
for i = 1:N
  if checkifcorrect 
      %all samples will be the same input, from UoB_ID
      samples(i,:) = CorrectCW;
  else 
     %these will be different random input binary sequences
     samples(i,:) = mod(input(i,:)*G,2);   
  end 
end  

%make steps starting from  0 to 0.2 inclusive, 30 intervals 
Prob_error = linspace(0,0.2,40);

%controls the error probability
for i = 1:length(Prob_error)
    %controls the number of simulations
    for j = 1:N
        
        %generate random errors using the error probability p
        %returns a sequence that may contain an error / bit flip
        REG(j,:)=+(rand(1,nout) < Prob_error(i));
        %add the error to the codeword with a the error
        codeword_errors(j,:) = mod(REG(j,:)+ samples(j,:),2);
        %recover using possible syndromes of codeword with errors
        s_det(j,:) = rem(codeword_errors(j,:)*H',2);
        for k = 1:8
            if (s_det(j,:)==syndrome(k,:))
                %check the table to look up which pattern it is
                %then find the 
                recovered_CW(j,:)=mod((table(k,:))+ codeword_errors(j,:),2);
            end
        end 
    end    
    %cacluation of Bit Error Rate - only looking at the first 4 bits which
    %are the most useful
    BER(i) = 1-length(find(recovered_CW(:,1:4) == Codeword(:,1:4)))/(4*N);
end

%theoretical BER
[PG,PDEC] = calculate_BER(nout,Prob_error);


%plots
semilogy(Prob_error,PG,'r','linewidth',1.5);
hold on
semilogy(Prob_error,BER,'b','linewidth',1.5);
hold on
semilogy(Prob_error,PDEC,'g','linewidth',1.5);
grid on;
legend('Theoretical BER with 1/2','Simulation BER','Theoretical BER with 3/15')
xlabel('Probability of Error')
ylabel('Error Rate')
title('Probability of Error & BER of (15,11) Hamming Code')




%% Functions and Utils

function [G,H] = hamming_15_11_G_H_matrices()

%G matrix has dimentions of n columns 
%G matrix has dimentions of k rows 
%Construct the Generator Matrix
G = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0; 
   	0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1];

%Take the dimentions of the generator matrix
[Row, Col] = size(G);
%Construction of the Parity matrix
%Take all rows and only n-k = 3 columns (take the parity matrix only)
%Take the part after the identity matrix
P = G(:,Row+1:Col);
%Transpose the parity matrix
P_T = P';
%Take the dimentions of the transposed Parity Matrix 
[Row_T,Col_T] = size(P_T);
%Construct an Identity of size of the row of the Transposed parity matrix
I = eye(Row_T);

fprintf("Task 1.1:\n")
%take the transpose of the parity matrix and pad the identity matrix
H = [P_T, I];

%1 function to reduce flooding of unnecessary variable



end

function Codeword = UOB_encoding_HC_UoBID(UoB_ID,G)

%convert it to decimal and then to a string
Binary_UoB_ID = dec2bin(UoB_ID);
%calculate lenght of the new decimal number - returns the number of characters in str.
lenghtofbinary = strlength(Binary_UoB_ID);
%take the 11 most right digits , str,start pos,end poss extractBetween(str,startPat,endPat)
encoder_input = extractBetween(Binary_UoB_ID,lenghtofbinary-10,lenghtofbinary);

%convert the cell from string to character array- put the 11 rightmost bits
encoder_input2 = char(encoder_input);
%get the codeword 
Codeword = mod(encoder_input2*G,2);

fprintf("Task 1.2:\n")
UoB_ID; %check

end

function [table,syndrome] = table_syndrome_construction(nout,kin,H)

%construct table showing you where the error is 
table = [zeros(1,15);eye(15)];


%15-11 = 4   => number of columns 
%n = 2^4 = 16  => number of rows 
    for n = 1:2^(nout-kin)
        %construct the syndrome table
        syndrome(n,:) = (table(n,:)*H');
    end
    
end

function [Codeword_out,Fixcodeword,CorrectCW] = decoding_codeword(Codeword,UoB_ID,H,table,syndrome)

%input codeword 
CorrectCW = Codeword;
%Adding Noise , flipping 1 bit IN THIS CASE: 
%stats counting the array from the left 
error_index = rem(UoB_ID,15)+1;
%save new codeword
Codeword_out = Codeword;
%apply the error in the codeword - output of the channel
Codeword_out(error_index) = ~Codeword_out(error_index);

%find the syndrome
s_1 = rem(Codeword_out*H',2);

%locate it in the table in the corresponding row
    for ns=1:8
        if (s_1 == syndrome(ns,:)) 
            index_ns = ns;
        end 
    end

%find from the table where the error is, output the error pattern 
error_pattern = (table(index_ns,:));
%find the position of the error
error_pos = find(error_pattern == 1);
%flip the error to find the correct codeword
Fixcodeword = Codeword_out;
Fixcodeword(error_pos) = ~Codeword_out(error_pos);

fixed_datainput = rem(Fixcodeword*H',2); %check if syndrome is 0000

end

function [PG,PDEC] = calculate_BER(nout,Prob_error)

    %using formula P(r<t)= Σ(n!/(n-r)!r!) x p^r (1-p)^n-r
    for r_errors=0:1
        total(r_errors+1,:)=(factorial(nout)/(factorial(nout-r_errors)*factorial(r_errors))).*Prob_error.^(r_errors).*(1-Prob_error).^(nout-r_errors);
    end 

%using formula 1-Σ of previous for general case assumption
PG=(1/2)*(1-(total(1,:)+total(2,:)));
%using formula 1-Σ of previous for Pdec Hamming Code
PDEC=(3/15)*(1-(total(1,:)+total(2,:)));

end

