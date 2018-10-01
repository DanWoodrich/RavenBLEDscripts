mooringsVec = dir('F:\Mid_Bering_mat_files\');

mooringsVec = mooringsVec(3:end);

for p = 1:size(mooringsVec,1)
    load(mooringsVec(p).name)
    R = load(mooringsVec(p).name);
    L = size(resltMTX);
    sumVecRW(p) = numel(find(resltMTX(:,:,1)==1));
    sumVecGS(p) = numel(find(resltMTX(:,:,3)==1));
    sumVecRWaGS(p) = numel(intersect(find(resltMTX(:,:,1)==1),find(resltMTX(:,:,3)==1)));
    B=0;
    C=0;
    D=0;
    E=0;
    for q = 2:L(3);
        B = vertcat(B,find(resltMTX(:,:,q)==1));
        C = vertcat(C,find(resltMTX(:,:,q)==2));
        BC = vertcat(B,C);
    end
        for r = [1:2 4:L(3)];
            D = vertcat(D,find(resltMTX(:,:,r)==1));
            E = vertcat(E,find(resltMTX(:,:,r)==2));
            DE = vertcat(D,E);
        end
    OnlyRW(p) = numel(setdiff(find(resltMTX(:,:,1)==1),BC));
    OnlyGS(p) = numel(setdiff(find(resltMTX(:,:,3)==1),DE));
end



mV2 = {mooringsVec(:).name};
finalVec = rot90([num2cell(sumVecRW); num2cell(sumVecGS); num2cell(sumVecRWaGS); num2cell(OnlyRW); num2cell(OnlyGS); mV2]);

T = cell2table(finalVec,'VariableNames',{'numRWy','numGSy','numRWGSy','numOnlyRW','numOnlyGS','Mooring'});

writetable(T,'F:\ALLRW_Counts.csv')