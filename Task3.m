clear all; close all; clc
%%TO SHOW ANSWERS REMOVE THE ; FROM THE END OF THE FUNCTIONS


%% Task 3: Convolutional Code

%Empirical investigation of the error correcting performance of a binary convolutional code

%number of samples for monte-carlo
N= 10000; 
k = 100;

%input is 100 
%2 outputs + redundant bits = 200 minimum

%Create the 4 states 
States = {   1   3     [0 0] [1 1]; ... state 00
             1   3     [1 1] [0 0]; ... state 01
             2   4     [0 1] [1 0]; ... state 10
             2   4     [1 0] [0 1]  ... state 11
    }; 

%Sorted as: next state if input is 0 followed by the output and 
%if input is 1 followed by the output

%               col 1:  next state if input = 0
%               col 2:  next state if input = 1
%               col 3:  output value if input = 0
%               col 4:  output value if input = 1

current_state = States(1,:);                % begins with state 1

%The σ is computed based on the required Eb/N0 or SNR (both can be used)
%We want to simulate what is the performance of the code at Eb/N0
%We first translate it to the linear scale 
%this will tell us with the assumption of signaling 
%what the signal noise should be  

%make steps starting from  0 to 0.2 inclusive, create an Eb/No array
ebn0db = linspace(0,10,21)



%adjust Eb/N0 for code rate, Increase  to simulate code rate for fair comparison
ebn0_adjusted = ebn0db + 10*log10(1/2);
%convert from db to normal value 
ebn0 = 10.^(ebn0db/10);
L_ebn0db = length(ebn0db)



for i= 1:L_ebn0db
    %inilitialise variables for calculates 
    error      = 0;
    u_error    = 0;
    u_e_adjust = 0;
    %calculate variance of ebn0
    sigma = 1 /(2*ebn0(i));
    %calculate standard deviation
    Sd = sqrt(sigma);
   
    for monte = 1:samples
        %generate random inputs
        input = genrate(k);
        
        %encode
        %output =
        
        
        
    end
    
end




%plot figures
%Task 1.1
figure(1);
title('Task 3');
xlabel('EbNo (db)');
ylabel('BER');
%semilogy();
hold on;



%% Functions and Utils

%generate awgn noise array
function noise = awgn(Sd, k, )
    noise(1,:) = Sd*randn(1,(k+3));
    noise(2,:) = Sd*randn(1,(k+3));
end

%generate erroneous codeword
function codeword = error(noise,output)
    for x = 1:2
        codeword(x,:) = output(x,:) + noise(x,:);
    end
end

%random input generation
function input = generate(k) 
    %generation of random input
    input = round(rand(1, k)); 
    %convert to BPSK (-1 and 1 instead of 0 and 1)
    input = input-1; 
    input = 2*input;
    %hardcode last 2 bits
    %terminating_bits=
    %input= [input, terminating_bits];
    
end



