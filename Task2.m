clear all; close all; clc
%%TO SHOW ANSWERS REMOVE THE ; FROM THE END OF THE FUNCTIONS

%%TASK2 Cyclic Redundancy Check (CRC) code

%% Task 2.0: CRC-5-USB code

%31 rows
nout = 31;
%26 columns 
kin = 26;
q = nout - kin; %31 - 26 = 5
%CRC-5-USB polynomial  
g=[1,0,0,1,0,1];
%all GFq polynomial operations in matlab 
%Therefore use reverse order mirror flip matrix g
gr = g(end:-1:1);
%gr = [1 0 1 0 0 1]

%checking if the generator polynomial is a valid (31,26) cyclic code
check_cyclic_code(gr);

%% Task 2.1: Syndrome caclulation

SyndromeUOB = CalcSyndromeUOB(nout,gr)



%% Task 2.2 
%Numerical instigation of error detecting capabilities of the CRC-5-USB code
%creates a random binary sequence as a matrix
%values from 0 to 1, rows 1, columns kin26 from 0 to 26
input = randi([0 1],1,kin);
%adds to the matrix input zeros at the end of the first row q=5 times.
d = [input,zeros(1,q)];
%mirror flip matrix d , the input matrix 
dr = d(end:-1:1);
%get the remainder and quot
[a,b] = gfdeconv(dr,gr,2);
%b = remainder
%flip the remainder
remainderf = b(end:-1:1);

%create a 1x5 array using extra zeros at the start + remainder at the end
check_value=[zeros(1,q-length(remainderf)),remainderf]; 
%add the 5bits at the end of the coded data
coded_data=[input,check_value];
number = 1:nout;

N = 1000;
number=1:nout;
error = 0;

%from 0 to 31
for hamming_weight=1:nout 
    %Monte Carlo simulation
    for j=1:N
        
        e = zeros(1,nout);% initial error sequence array of 1x31 zeros
        rand_index = randperm(nout); % random rank the index
        draw_rand_index = rand_index(1:hamming_weight); % take the m numbers out
        flip = number(draw_rand_index);
        e(flip) = 1;
        encoded_data=mod((coded_data+e),2);
        E_dataf = encoded_data(end:-1:1);
        [a1,b1]=gfdeconv(E_dataf,gr,2);
        %a1r = a1(end:-1:1);
        b1r = b1(end:-1:1);
        b1_stuff=[zeros(1,(q-length(b1r))),b1r];
        if (length(find(b1_stuff==0))==5) % judge whether the error is detetable
            error=error+1; %increment error
        end
    end
    X(hamming_weight)=hamming_weight; %calcualting hamming weight of the error in an array
    P(hamming_weight)=1-error/N; %calculating detection probability in an array
    error=0; %resetting the error
end

bar(X,P) %plot in a bar graph
title('Detection Capability of the CRC-5-USB (31,26) code as a function of Hamming weight of the error')
xlabel('Hamming weight of the error')
ylabel('Detection probaility')




%% Functions and Utils

function syndrome = CalcSyndromeUOB(nout,gr)
    
    %Input = uob number 1637541 / 0805547 / 1145506
    uob_id = 1637541;
    uob_bin = dec2bin(uob_id);
    lenghtofbinary = strlength(uob_bin); %length of string
    %convert from string to array  
    for i=1: lenghtofbinary %loop over the string's length
        uob_array(i) = [str2double(uob_bin(i))]; %converter
    end

    %add zeros infront 
    uob_zeros = [zeros(1,nout-lenghtofbinary), uob_array];
    length_uob_zeros = length(uob_zeros); %check length of array
    %flipping array 
    input_array_r = uob_zeros(end:-1:1);
    
    %Encoding is done as as simple polynomial multiplication c(x)=a(x)g(x)
    %Therefore for deconvolution the reverse process had to be done
    %To find syndrome take the remainder of the deconvolution with g matrix
    [a,b] = gfdeconv(input_array_r,gr,2);
    syndrome = b
    syndrome = b(end:-1:1);
    
    fprintf("Task 2.1: ")

end


function remainder = check_cyclic_code(gr)

%checking if g(x) is a polynomial of a valid cyclic code
%i.e. is g(x) a factor of p(x) = x^n+1  p(x) = (x^31)+1
p =[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
%flip i t
pr = p(end:-1:1);
%so if we take the polynomial p(x) by g(x) and there is no remainder after
%the division then that g(x) is a valid generator matrix
%The multiplication of polynomials is a obtained by the convolution
%Division of polynomials is obtained by the deconvolution
[quot,remd] = gfdeconv(pr,gr);
%since remainder = 0 this generator polynomial is a valid (31,26) cyclic code
%In this case
%quot [1,0,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1] 1x27
%remd 0 

fprintf("Task 2.0: \n")
if remd == 0
    fprintf("Valid \n\n")
else
    fprintf("Invalid \n\n")

end

end



