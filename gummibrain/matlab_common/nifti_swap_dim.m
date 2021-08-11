function nifti_swap_dim(infile,outfile,neworder)



if(numel(neworder) == 3)
    neworder(4)=4;
end

M=load_nifti(infile);
M.dim([2 3 4 5])=M.dim(neworder+1);
M.vol=permute(M.vol,neworder);

save_nifti(M,outfile);
