clear all; close all; clc;


%% Ideal case
clear all; close all; clc; 
load('Xt1_1.mat');
load('Xt2_1.mat');
load('Xt3_1.mat');

% Each of these pairs of vectors will be different lengths, so what should
% you do to not get dimension errors?

Xt2_1 = Xt2_1(:,1:226);
Xt3_1 = Xt3_1(:,1:226);
%figure(3)
%surf(Xt1_1,Xt2_1,Xt3_1)

% After taking care of the dimension issue, subtract the mean off of each
% component; i.e. each row vector for each video (see documentation for the
% "mean" function and becareful of the dimensions to make sure it's done
% correctly because it's finicky).

means1 = mean(Xt1_1(1,:));
means2 = mean(Xt1_1(2,:));
Xt1_1(1,:) = Xt1_1(1,:) - repmat(means1,1,226);
Xt1_1(2,:) = Xt1_1(2,:) - repmat(means2,1,226);

means1 = mean(Xt2_1(1,:));
means2 = mean(Xt2_1(2,:));
Xt2_1(1,:) = Xt2_1(1,:) - repmat(means1,1,226);
Xt2_1(2,:) = Xt2_1(2,:) - repmat(means2,1,226);


means1 = mean(Xt3_1(1,:));
means2 = mean(Xt3_1(2,:));
Xt3_1(1,:) = Xt3_1(1,:) - repmat(means1,1,226);
Xt3_1(2,:) = Xt3_1(2,:) - repmat(means2,1,226);


% Combine the pairs of vectors into one 6xL matrix, where L is the length,
% which will vary depending on the video.
% Take the svd and do a coordinate transformation to the principal
% components just like we did in the lectures (we called this matrix Y).
% Save the matrix Y as A1.

CombMtrx = [Xt1_1; Xt2_1; Xt3_1];
[U,S,V] = svd(CombMtrx);
Y = U*CombMtrx;

A1 = Y; % Shape: 6x226 double


% Save the energies of nontrivial the singular values as a vector in the
% varialbe A2.
A2 = zeros(6,1); % Shape: 6x1 double
sig = diag(S);
A2(1,1) = sig(1)^2/sum(sig.^2);
A2(2,1) = sum(sig(1:2).^2)/sum(sig.^2);
A2(3,1) = sum(sig(1:3).^2)/sum(sig.^2);
A2(4,1) = sum(sig(1:4).^2)/sum(sig.^2);
A2(5,1) = sum(sig(1:5).^2)/sum(sig.^2);
A2(6,1) = sum(sig(1:5).^2)/sum(sig.^2);

%%% Plot the energies with a log-scaled ordinate for the report (not for 
% the autograder)

%figure(1)
%semilogy(sig,'o','LineWidth',2)
%title('Test 1 Singular Value Energies on Log-Scaled Ordinate')
%ylabel('Energy')
%xlabel('Singular Value')
%grid on

% The appropriate rank necessary for a decent approximation will be
% informed by the energies and by judging the qualitative convergence (no
% need for rigorous analysis, you can just plot and check) as we increase
% the rank of the approximation.
% Save the rank-n approximation as A3.
% Of course, this will be subjective (for now), so you may have to try it
% a few times to match what I picked.  You also may not agree with my
% choice, or your mass tracking may be different based on what you clicked.
% If you get it to pass the autograder, but you disagree put it in the
% report.
A3 = U(:,1)*S(1,1)*V(:,1)' + U(:,2)*S(2,2)*V(:,2)'; % Shape: 6x226 double
%plot3(A3(1,:),A3(2,:),A3(3,:),'g.','Markersize',10)
grid on;

%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence

x = linspace(0,1,6);
t = linspace(0,2,226);
[T, X] = meshgrid(t,x);

%figure(2)
%surf(X,T,A3)
%title('Rank 2 Approximation of Test 1')
%shading interp



%% Test 2
clear all; close all; clc;
load('Xt1_2.mat')
load('Xt2_2.mat')
load('Xt3_2.mat')

% Repeat what you did for the Ideal Case
% Each of these pairs of vectors will be different lengths, so what should
% you do to not get dimension errors?

Xt2_2 = Xt2_2(:,1:314);
Xt3_2 = Xt3_2(:,1:314);




% After taking care of the dimension issue, subtract the mean off of each
% component; i.e. each row vector for each video (see documentation for the
% "mean" function and becareful of the dimensions to make sure it's done
% correctly because it's finicky).

means1 = mean(Xt1_2(1,:));
means2 = mean(Xt1_2(2,:));
Xt1_2(1,:) = Xt1_2(1,:) - repmat(means1,1,314);
Xt1_2(2,:) = Xt1_2(2,:) - repmat(means2,1,314);

means1 = mean(Xt2_2(1,:));
means2 = mean(Xt2_2(2,:));
Xt2_2(1,:) = Xt2_2(1,:) - repmat(means1,1,314);
Xt2_2(2,:) = Xt2_2(2,:) - repmat(means2,1,314);


means1 = mean(Xt3_2(1,:));
means2 = mean(Xt3_2(2,:));
Xt3_2(1,:) = Xt3_2(1,:) - repmat(means1,1,314);
Xt3_2(2,:) = Xt3_2(2,:) - repmat(means2,1,314);



% Combine the pairs of vectors into one 6xL matrix, where L is the length,
% which will vary depending on the video.
% Take the svd and do a coordinate transformation to the principal
% components just like we did in the lectures (we called this matrix Y).
% Save the matrix Y as A4.

CombMtrx = [Xt1_2; Xt2_2; Xt3_2];
[U,S,V] = svd(CombMtrx);
Y = U*CombMtrx;

A4 = Y; % 6x314 double


% Save the energies of nontrivial the singular values as a vector in the
% varialbe A5.
A5 = zeros(6,1); % Shape: 6x1 double
sig = diag(S);
A5(1,1) = sig(1)^2/sum(sig.^2);
A5(2,1) = sum(sig(1:2).^2)/sum(sig.^2);
A5(3,1) = sum(sig(1:3).^2)/sum(sig.^2);
A5(4,1) = sum(sig(1:4).^2)/sum(sig.^2);
A5(5,1) = sum(sig(1:5).^2)/sum(sig.^2);
A5(6,1) = sum(sig(1:6).^2)/sum(sig.^2);


%%% Plot the energies with a log-scaled ordinate for the report (not for 
% the autograder)



% The appropriate rank necessary for a decent approximation will be
% informed by the energies and by judging the qualitative convergence (no
% need for rigorous analysis, you can just plot and check) as we increase
% the rank of the approximation.
% Save the rank-n approximation as A6.
% Of course, this will be subjective (for now), so you may have to try it
% a few times to match what I picked.  You also may not agree with my
% choice, or your mass tracking may be different based on what you clicked.
% If you get it to pass the autograder, but you disagree put it in the
% report.
A6 = U(:,1)*S(1,1)*V(:,1)'+ U(:,2)*S(2,2)*V(:,2)' + U(:,3)*S(3,3)*V(:,3)';% + U(:,4)*S(4,4)*V(:,4)'+ U(:,5)*S(5,5)*V(:,5)'; % 6x314 double


%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence

%figure(1)
%semilogy(sig,'o','LineWidth',2)
%title('Test 2 Singular Value Energies on Log-Scaled Ordinate')
%ylabel('Energy')
%xlabel('Singular Value')
%grid on
%
%
%x = linspace(0,1,6);
%t = linspace(0,2,314);
%[T, X] = meshgrid(t,x);
%
%figure(2)
%surf(X,T,A6)
%title('Rank 3 Approximation of Test 2')
%shading interp




%% Test 3

clear all; close all; clc;
load('Xt1_3.mat')
load('Xt2_3.mat')
load('Xt3_3.mat')

% Repeat what you did for the Ideal Case
% Each of these pairs of vectors will be different lengths, so what should
% you do to not get dimension errors?
Xt1_3 = Xt1_3(:,1:237);
Xt2_3 = Xt2_3(:,1:237);





% After taking care of the dimension issue, subtract the mean off of each
% component; i.e. each row vector for each video (see documentation for the
% "mean" function and becareful of the dimensions to make sure it's done
% correctly because it's finicky).


means1 = mean(Xt1_3(1,:));
means2 = mean(Xt1_3(2,:));
Xt1_3(1,:) = Xt1_3(1,:) - repmat(means1,1,237);
Xt1_3(2,:) = Xt1_3(2,:) - repmat(means2,1,237);

means1 = mean(Xt2_3(1,:));
means2 = mean(Xt2_3(2,:));
Xt2_3(1,:) = Xt2_3(1,:) - repmat(means1,1,237);
Xt2_3(2,:) = Xt2_3(2,:) - repmat(means2,1,237);


means1 = mean(Xt3_3(1,:));
means2 = mean(Xt3_3(2,:));
Xt3_3(1,:) = Xt3_3(1,:) - repmat(means1,1,237);
Xt3_3(2,:) = Xt3_3(2,:) - repmat(means2,1,237);




% Combine the pairs of vectors into one 6xL matrix, where L is the length,
% which will vary depending on the video.
% Take the svd and do a coordinate transformation to the principal
% components just like we did in the lectures (we called this matrix Y).
% Save the matrix Y as A7.
CombMtrx = [Xt1_3;Xt2_3;Xt3_3];
[U,S,V] = svd(CombMtrx);
Y = U*CombMtrx;
A7 = Y; % Shape: 6x237 double

sig = diag(S);

% Save the energies of nontrivial the singular values as a vector in the
% varialbe A8.
A8 = zeros(6,1);% Shape: 6x1 double
A8(1,1) = sig(1)^2/sum(sig.^2);
A8(2,1) = sum(sig(1:2).^2)/sum(sig.^2);
A8(3,1) = sum(sig(1:3).^2)/sum(sig.^2);
A8(4,1) = sum(sig(1:4).^2)/sum(sig.^2);
A8(5,1) = sum(sig(1:5).^2)/sum(sig.^2);
A8(6,1) = sum(sig(1:6).^2)/sum(sig.^2);


%%% Plot the energies with a log-scaled ordinate for the report (not for 
% the autograder)






% The appropriate rank necessary for a decent approximation will be
% informed by the energies and by judging the qualitative convergence (no
% need for rigorous analysis, you can just plot and check) as we increase
% the rank of the approximation.
% Save the rank-n approximation as A9.
% Of course, this will be subjective (for now), so you may have to try it
% a few times to match what I picked.  You also may not agree with my
% choice, or your mass tracking may be different based on what you clicked.
% If you get it to pass the autograder, but you disagree put it in the
% report.
A9 = U(:,1)*S(1,1)*V(:,1)' + U(:,2)*S(2,2)*V(:,2)' + U(:,3)*S(3,3)*V(:,3)'; % 6x237 double

%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence

%figure(1)
%semilogy(sig,'o','LineWidth',2)
%title('Test 3 Singular Value Energies on Log-Scaled Ordinate')
%ylabel('Energy')
%xlabel('Singular Value')
%grid on
%
%
%x = linspace(0,1,6);
%t = linspace(0,2,237);
%[T, X] = meshgrid(t,x);
%
%figure(2)
%surf(X,T,A9)
%title('Rank 3 Approximation of Test 3')
%shading interp






%% Test 4
clear all; close all; clc;
load('Xt1_4.mat')
load('Xt2_4.mat')
load('Xt3_4.mat')
%
% Repeat what you did for the Ideal Case
% Each of these pairs of vectors will be different lengths, so what should
% you do to not get dimension errors?

Xt2_4 = Xt2_4(:,1:392);
Xt3_4 = Xt3_4(:,1:392);


% After taking care of the dimension issue, subtract the mean off of each
% component; i.e. each row vector for each video (see documentation for the
% "mean" function and becareful of the dimensions to make sure it's done
% correctly because it's finicky).


% Combine the pairs of vectors into one 6xL matrix, where L is the length,
% which will vary depending on the video.
% Take the svd and do a coordinate transformation to the principal
% components just like we did in the lectures (we called this matrix Y).
% Save the matrix Y as A10.
CombMtrx = [Xt1_4;Xt2_4;Xt3_4];
[U,S,V] = svd(CombMtrx);
Y = U*CombMtrx;
A10 = Y; % 6x392 double


% Save the energies of nontrivial the singular values as a vector in the
% varialbe A11.
sig = diag(S);
A11 = zeros(6,1);% 6x1 double

A11(1,1) = sig(1)^2/sum(sig.^2);
A11(2,1) = sum(sig(1:2).^2)/sum(sig.^2);
A11(3,1) = sum(sig(1:3).^2)/sum(sig.^2);
A11(4,1) = sum(sig(1:4).^2)/sum(sig.^2);
A11(5,1) = sum(sig(1:5).^2)/sum(sig.^2);
A11(6,1) = sum(sig(1:6).^2)/sum(sig.^2);

%%% Plot the energies with a log-scaled ordinate for the report (not for 
% the autograder)






% The appropriate rank necessary for a decent approximation will be
% informed by the energies and by judging the qualitative convergence (no
% need for rigorous analysis, you can just plot and check) as we increase
% the rank of the approximation.
% Save the rank-n approximation as A12.
% Of course, this will be subjective (for now), so you may have to try it
% a few times to match what I picked.  You also may not agree with my
% choice, or your mass tracking may be different based on what you clicked.
% If you get it to pass the autograder, but you disagree put it in the
% report.
A12 = U(:,1)*S(1,1)*V(:,1)' + U(:,2)*S(2,2)*V(:,2)' ;% U(:,3)*S(3,3)*V(:,3)';% + U(:,4)*S(4,4)*V(:,4)'; % 6x392 double


%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence



figure(1)
semilogy(sig,'o','LineWidth',2)
title('Test 4 Singular Value Energies on Log-Scaled Ordinate')
ylabel('Energy')
xlabel('Singular Value')
grid on


x = linspace(0,1,6);
t = linspace(0,2,392);
[T, X] = meshgrid(t,x);

figure(2)
surf(X,T,A12)
title('Rank 4 Approximation of Test 4')
shading interp

