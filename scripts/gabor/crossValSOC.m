function [spaceTrain,allTrain,allTest] = crossValSOC()
%
% function generates stimuli for 10 fold cross validalidation, with the
% space/zebra stimuli appearing in the training set.
%
% Example:
% [spaceTrain,allTrain,allTest] = crossValSOC()

% These are the space & zebra stimuli:
stimSpace = [1:38 59:68 74:78];
% Devide in 10 for 10-fold crossvalidation, first randomize:
randSpace = stimSpace(randperm(length(stimSpace)));
% test                              % train
spaceSet{1,1} = randSpace(1:5);     spaceSet{1,2} = randSpace(6:end);
spaceSet{2,1} = randSpace(6:10);    spaceSet{2,2} = randSpace([1:5 11:end]);
spaceSet{3,1} = randSpace(11:15);   spaceSet{3,2} = randSpace([1:10 16:end]);
spaceSet{4,1} = randSpace(16:20);   spaceSet{4,2} = randSpace([1:15 21:end]);
spaceSet{5,1} = randSpace(21:25);   spaceSet{5,2} = randSpace([1:20 26:end]);
spaceSet{6,1} = randSpace(26:30);   spaceSet{6,2} = randSpace([1:25 31:end]);
spaceSet{7,1} = randSpace(31:35);   spaceSet{7,2} = randSpace([1:30 36:end]);
spaceSet{8,1} = randSpace(36:41);   spaceSet{8,2} = randSpace([1:35 42:end]);
spaceSet{9,1} = randSpace(42:47);   spaceSet{9,2} = randSpace([1:41 48:end]);
spaceSet{10,1} = randSpace(48:53);  spaceSet{10,2} = randSpace([1:47]);

% These are the non-space stimuli:
stimOther = [39:58 69:73 79:86];
% Devide in 10 for 10-fold crossvalidation, first randomize:
randOther = stimOther(randperm(length(stimOther)));
% test                              % train
otherSet{1,1} = randOther(1:3);     otherSet{1,2} = randOther(4:end);
otherSet{2,1} = randOther(4:6);     otherSet{2,2} = randOther([1:3 7:end]);
otherSet{3,1} = randOther(7:9);     otherSet{3,2} = randOther([1:6 10:end]);
otherSet{4,1} = randOther(10:12);   otherSet{4,2} = randOther([1:9 13:end]);   
otherSet{5,1} = randOther(13:15);   otherSet{5,2} = randOther([1:12 16:end]);
otherSet{6,1} = randOther(16:18);   otherSet{6,2} = randOther([1:15 19:end]);
otherSet{7,1} = randOther(19:21);   otherSet{7,2} = randOther([1:18 22:end]);
otherSet{8,1} = randOther(22:25);   otherSet{8,2} = randOther([1:21 26:end]);
otherSet{9,1} = randOther(26:29);   otherSet{9,2} = randOther([1:25 30:end]);
otherSet{10,1} = randOther(30:33);  otherSet{10,2} = randOther([1:29]);

spaceTrain = {};
allTrain = {};
allTest = {};
for kk = 1:10
    spaceTrain{kk} = sort(spaceSet{kk,2});
    allTrain{kk} = sort([spaceSet{kk,2} otherSet{kk,2}]);
    allTest{kk} = sort([spaceSet{kk,1} otherSet{kk,1}]);
end


