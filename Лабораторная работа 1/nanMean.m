function meanXMat = nanMean(XMat)
    aMat = ones(size(XMat));
    aMat(isnan(XMat)) = 0;
    XMat(isnan(XMat)) = 0;
    meanXMat = sum(XMat)./sum(aMat);
end