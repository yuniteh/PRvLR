function feat = addNoise(daqdata, ratio)
frame = 199;
inc = 25;
feat_val = 15;
feat_num = 4;
num_emg = size(daqdata, 2);
num_samps = size(daqdata, 1);
s_down = 1:inc:num_samps - frame;           % incremented window indices
feat = zeros(length(s_down),num_emg*feat_num);
noise = ratio * randn(num_samps, 1);
noisy_data = zeros(num_samps, num_emg);
figure 
for i = 1:num_emg
    subplot(num_emg,1,i)
    plot(daqdata(:,i))
end
daqdata = (double(daqdata)/((2^16-1)/10)) - 5;

figure 
for i = 1:num_emg
    subplot(num_emg,1,i)
    plot(daqdata(:,i))
end

figure
for i = 1:num_emg
    noisy_data(:,i) = daqdata(:,i) + noise;
    subplot(num_emg,1,i)
    plot(noisy_data(:,i))
    ylim([-5,5])
end
% convert to uint16
noisy_data = uint16((noisy_data + 5) * (2^16 - 1)/10);

for i = 1:length(s_down)
    feat(i,:) = TDFeatExtractmex(feat_val,noisy_data(s_down(i):s_down(i)+frame,:)');
end

end