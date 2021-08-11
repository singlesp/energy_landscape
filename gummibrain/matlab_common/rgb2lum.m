function L = rgb2lum(rgb, displaycalib)

lr = displaycalib.aR*rgb(:,1).^( displaycalib.bR);
lg = displaycalib.aG*rgb(:,2).^( displaycalib.bG);
lb = displaycalib.aB*rgb(:,3).^( displaycalib.bB);

L = lr+lg+lb+displaycalib.L0;