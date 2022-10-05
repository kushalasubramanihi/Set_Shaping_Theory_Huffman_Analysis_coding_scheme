%              Set_Shaping_Theory_Huffman_Analysis_coding_scheme
%
%              email:kushala.subramanihi@outlook.com
%
% In this version of the program presented in the following article:
%
% "Practical applications of Set Shaping Theory in Huffman coding"
% https://arxiv.org/abs/2208.13020
%
% I introduce an analysis about the number of bits needed to reproduce 
% Huffman tree (decoding scheme). This data is fundamental, because it serves
% to understand that the gain obtained by applying the set shaping theory 
% does not imply an increase in the number of bits necessary to represent 
% the decoding scheme (correspondence between codeword symbols).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Definition
%
% Definition: given a sequence S, we call Cs(S) the coded sequence in which 
% the symbols are replaced by a uniquely decodable code. 
%
% Definition: given a sequence S of random variables of length N over the 
% alphabet A=(s1,s2.....sh), we call the zero order empirical entropy H0(S) 
% the function defined by:
% 
%  H0(S)=-∑p(si)log2p(si)=-p(s1)log2p(s1).................-p(sh)log2p(sh)
%  
%  with p(si) denote the number of occurrences of the symbol si inside S.
%
% The zero order empirical entropy  is a lower bound for a compressor which
% encodes each symbol independently of any of its preceding or following symbols.
% Therefore, NH0(S) represents the output size of an ideal compressor.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Analysis of the coding scheme (Huffman tree)
%
% The first calculated parameter is defined by the sum of all the codewords.
% This parameter gives us an indicative value of the number of bits needed 
% to represent the dictionary (list of codewords) through a simple table 
% of correspondence codeword symbols.
%
% Using a correspondence table is not an effecient way to represent the 
% condig scheme. Much research has been done to find the most efficient way 
% to represent the coding scheme. In the following article "T. Gagie, G. Navarro, 
% Y. Nekrich, and A. O. Pereira. Efficient and compact representations of 
% prefix codes. IEEE Trans. Inf. Theory, 61 (9): 4999–5011,2015." the author 
% gave a representation of a coding capable to encode or decode a symbol 
% in constant worst case time. This code uses |A|lglmax + o(|A|) + O(lmax^2) 
% bits of space, where |A| and lmax are the alphabet size and maximum 
% codeword length. Therefore, the parameters |A| and lmax turn out 
% to be a good indicator of how much the coding scheme can be compressed.
%
% In conclusion these three parameters are calculated both for the generated 
% sequence and for the transformed sequence:
%
% 1) sum of all the codewords
% 2) lmax maximum codeword length
% 3) |A| alphabet size
%
% The results obtained show that the parameter relative to the sum of the 
% codewords and the parameter of the dimension of the alphabet decrease by 
% applying the set shaping theory, instead the parameter lmax remains constant.
%
% This results are very clear applying the set shaping theory does 
% not imply an increase in the number of bits necessary to describe the 
% coding scheme, on the contrary a slight decrease is observed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The program performs the following operations:
% 1) generates a random sequence S with uniform distribution
% 2) calculate the zero order empirical entropy H0(S)
% 3) apply the Huffman coding on S and calculate a series of fundamental parameters 
%    to understand the limit of compression of the coding scheme
% 4) the set shaping theory is performed, the sequence S is transformed into 
%    the sequence f(S)
% 5) apply the Huffman coding on f(S) and calculate a series of fundamental 
%    parameters to understand the limit of compression of the coding scheme
% 6) compares NH0(S) of the generated sequence with the length of the
%    encoded transformated sequence cs(f(S)) 
% 7) repeats all these steps a number of times defined by the parameter history
% 8) display the average values obtained
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Important
%
% If you change the length of the sequence and the number of the generated 
% symbols, you have to be careful that the Huffman encoding approximates the 
% coding limit of the sequence by about one bit. if you take too 
% long sequences the Huffman algorithm becomes very inefficient therefore, 
% it cannot detect the advantage obtained by the transformation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;

% Imput of the random sequence S

% As a general rule, if you take ns symbols the length of the sequence 
% must be about 2*ns, in this case the Huffman encoding approximates 
% the coding limit of the sequence by about one bit.

% Number of sequences generated
history=1000; 
% Number of symbols
ns=40; 
% Lenght of the sequence
len=80;  

% ns must be in this range 20 <= ns <= 500 
% len must be in this range 40 <= len <= 1000

cs=0;
totcodel=0;
totlc=0;
tottlc=0;
itnent=0;
lc=0;
tlc=0;
itlc=0;
itlensd=0;
ttlensd=0;
timaxcodeword=0;
ttmaxcodeword=0;
tima=0;
ttma=0;

for i=1:history

 % Generation of the sequence S with a uniform distribution

 symbols=1:ns;
 prob(1,1:ns)=1/ns;
 seq=randsrc(1,len,[symbols; prob]);

 % We calculate NH0(S) (coding limit of S)

 lc=0;

 for i2=1:len

  sy=seq(1,i2);
  fs=nnz(seq==sy)/len;
  lc=lc-log2(fs);

 end

 % we define the length of the vectors that are used in the encoding

 index=0;

 for i2=1:ns

  fs=nnz(seq==i2)/len;

  if fs > 0

   index=index+1;    

  end

 end 

 c=zeros(index,1);
 vs=zeros(index,1);
 index=0;

 % We calculate the frequencies of the symbols in the sequence S

 for i2=1:ns

  fs=nnz(seq==i2)/len;

  if fs > 0

   index=index+1;    
   c(index)=nnz(seq==i2)/len;  
   vs(index)=i2;

  end

 end 

 % We code the sequence S

 counts=c;

 idict=huffmandict(vs,counts);
 icode=huffmanenco(seq,idict);

 bcode=de2bi(icode);
 icodel=numel(bcode);

 % Analysis of the coding scheme (Hufmann coding applied to the initial sequence S)

 ilensd=0;
 imaxcodeword=0;
 sdict=size(idict);
 ima=sdict(1,1);

 for i2=1:index

  d=idict(i2,2);
  sd=numel(d{1});
  ilensd=ilensd+sd;

  if sd > imaxcodeword 

   imaxcodeword=sd;

  end

 end

 % Start trasformation f(S)

 mcodel=10000;

 nseq=fSSTt(seq);

 % The new sequence is long nlen=len+1
 
 nlen=len+1;

 % We calculate N2H0(f(S)) with N2=N+1 (coding limit of f(S))

  tlc=0;

  for i2=1:nlen

   sy=nseq(1,i2);
   fs=nnz(nseq==sy)/nlen;
   tlc=tlc-log2(fs);

  end 

 % Having transformed the sequence, we have to redefine the length of the 
 % vectors that are used in the encoding

 index=0;

 for i2=1:ns

  fs=nnz(nseq==i2)/nlen;

  if fs > 0

   index=index+1;    

  end

 end 

 c=zeros(index,1);
 vs=zeros(index,1);
 index=0;

 % We calculate the frequencies of the symbols in the transformed sequence f(S)

 for i2=1:ns

  fs=nnz(nseq==i2)/nlen;

  if fs > 0

   index=index+1;    
   c(index)=nnz(nseq==i2)/nlen;  
   vs(index)=i2;

  end

 end 

 % We code the transformed sequence f(S)

 counts=c;

 tdict=huffmandict(vs,counts);
 tcode=huffmanenco(nseq,tdict);

 bcode=de2bi(tcode);
 tcodel=numel(bcode);

 % Analysis of the coding scheme (Hufmann coding applied to the transformed sequence f(S))

 tlensd=0;
 tmaxcodeword=0;
 sdict=size(tdict);
 tma=sdict(1,1);

 for i2=1:index

  d=tdict(i2,2);
  sd=numel(d{1});
  tlensd=tlensd+sd;

  if sd > tmaxcodeword 

   tmaxcodeword=sd;

  end

 end

 % If the length of the encoded message of the transformed sequence is less
 % than coding limit NH0(S) of the original sequence S, we increase the counter
 % cs by one

 if tcodel < lc

  cs=cs+1;

 end

 % We apply the inverse transform and we obtain the initial sequence S

 iseq=invfSSTt(nseq);

 % We check that the obtained sequence is equal to the initial sequence S

 flag=isequal(seq,iseq);

 if flag == false

   fprintf('Error, sequence not equal to the initial sequence\n');

 end

 totcodel=totcodel+tcodel;
 totlc=totlc+lc;
 tottlc=tottlc+tlc;
 itlensd=itlensd+ilensd;
 ttlensd=ttlensd+tlensd;
 timaxcodeword=timaxcodeword+imaxcodeword;
 ttmaxcodeword=ttmaxcodeword+tmaxcodeword;
 tima=tima+ima;
 ttma=ttma+tma;

end

% we calculate the average values of all parameters
  
 medlc=totlc/history;
 medcodel=totcodel/history;
 medtlc=tottlc/history;
 mitlensd=itlensd/history;
 mttlensd=ttlensd/history;
 mimaxcodeword=timaxcodeword/history;
 mtmaxcodeword=ttmaxcodeword/history;
 mima=tima/history;
 mtma=ttma/history;

% We calculate the percentage of sequences where the length of the encoded
% transformed sequence Cs(f(S)) is less than the coding limit NH0(S) of the 
% generated sequence S

 pcs=(cs/history)*100;

% We display the average values obtained

 fprintf('The average of the coding limit NH0(S) of the generated sequences S\n');
 medlc

 fprintf('The average of the length of the encoded transformed sequence Cs(f(S)) \n');
 medcodel

 fprintf('The average of the coding limit N2H0(f(S)) of the transformed sequences f(S) with N2=N+1 \n');
 medtlc

 fprintf('The average lenght of the dictionary (sum of all the codewords, Hufmann coding applied to the initial sequence S)\n');
 mitlensd

 fprintf('The average lenght of the dictionary (sum of all the codewords, Hufmann coding applied to the transformed sequence f(S))\n');
 mttlensd

 fprintf('The average maximum codeword length of the dictionary (Hufmann coding applied to the initial sequence S)\n');
 mimaxcodeword

 fprintf('The average maximum codeword length of the dictionary (Hufmann coding applied to the transformed sequence f(S)) \n');
 mtmaxcodeword

 fprintf('the average alphabet size initial sequence S\n');
 mima

 fprintf('the average alphabet size transformed sequence f(S)\n');
 mtma

 fprintf('Number of sequences where the length of the encoded transformed sequence Cs(f(S)) is less than the coding limit NH0(S) of the generated sequence S\n');
 cs

 fprintf('There is a percentage of %2.0f%% that length of the encoded transformed sequence Cs(f(S)) is less than the coding limit NH0(S) of the generated sequence S (N lenght of the sequence, H0(S) zero order empirical entropy)\n',pcs);

