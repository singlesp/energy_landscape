function blockperiods = BlockPeriods(v,minperiod)
 
S = (v==0);
d = diff([0; S; 0]);
s_start = find(d==1);
s_end = find(d==-1)-1;
s_dur = diff([s_start s_end],[],2);

silence_period = [s_start s_end s_dur];
silence_period = silence_period(silence_period(:,3)>minperiod,:);


M = ones(size(S));
for i = 1:size(silence_period,1)
    sp = silence_period(i,:);
    M(sp(1):sp(2)) = 0;
end
d = diff([0; M; 0]);
m_start = find(d==1);
m_end = find(d==-1)-1;
m_dur = diff([m_start m_end],[],2);

movie_period = [m_start m_end m_dur];

full_period = [[zeros(size(silence_period,1),1) silence_period]; 
    [(1:size(movie_period,1))' movie_period]];

blockperiods = sortrows(full_period,2);
