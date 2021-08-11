function Ymat = applyTEO(datamat, order, fs, lowcutoff, highcutoff)

Ymat = zeros(size(datamat));
for i = 1:size(datamat,2)
	data = datamat(:,i);
	if(nargin >= 3)
		data = amri_sig_filtfft(data,fs,lowcutoff,highcutoff); 
	end

	k = order;

	Ekm = ones(size(data))*nan;
	Ekp = ones(size(data))*nan;
	Ekm(k+1:end) = data(1:end-k);
	Ekp(1:end-k) = data(k+1:end);
	Y = data.^2-Ekm.*Ekp;
	Y(isnan(Y))=0;
	Y(Y<=0)=0;
	Ymat(:,i) = Y;
end