% Загружаем исходный файл

dir_name = sprintf('E:\Netlist\ITGLOBAL\wi_fi\matlab\');
file_name = sprintf('var.mat',dir_name);

var = importdata(file_name); 

%%
%  for i = 1:7120-16
%      a(i) =  waveStruct.waveform(i) * conj(waveStruct.waveform(i+16));
%  end
% 
%     p_i = waveStruct.waveform .* conj(waveStruct.waveform);
% 
%  for i = 1:106
%      a_sum(i) =  real(sum(a((i*16)-15:(i*16))));
%      p_sum(i) =  sum(p_i((i*16)-15:(i*16)));
%  end
% 
% r_sum = p_sum.';
% 
% corr = abs(a_sum) ./ r_sum;

%%
% rxWaveform = [randi([0, 1], 100, 1); var.waveform];
rxWaveform = [var.waveform];
offset = 0;
threshold = 0.75;

startOffset = wlanPacketDetect(rxWaveform, var.config.waveform.ChannelBandwidth, offset, threshold);
% totalOffset = offset + startOffset;

%% Числитель
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fifo 
delay_sample_fifo_out = [zeros(17,1); rxWaveform(2:end-16)];
% complex_mult_delay
delay_prod_inst = rxWaveform .* conj(delay_sample_fifo_out);

% avg_dual_channel
delay_prod = delay_prod_inst;

for i = 2 : size(delay_prod_inst)
    if (i < 34)
        delay_prod_inst(i) = delay_prod_inst(i) + delay_prod_inst(i-1);
    else
        delay_prod_inst(i) = delay_prod_inst(i-1) + delay_prod(i) - delay_prod(i-16);
    end
end

% abs
i_delay = abs(real(delay_prod_inst));
q_delay = abs(imag(delay_prod_inst));

% find max
for i = 1:size(delay_prod_inst)
    if (i_delay(i) > q_delay(i))
        max_d(i) = i_delay(i);
    else
        max_d(i) = q_delay(i);
    end
end

% find min
for i = 1:size(delay_prod_inst)
    if (i_delay(i) > q_delay(i))
        min_d(i) = q_delay(i);
    else
        min_d(i) = i_delay(i);
    end
end

% mag
alpha = 1;
beta = 1/4;

mag = alpha*max_d + beta*min_d;
mag = mag.';

% mag traditional
mag_double = sqrt(real(delay_prod_inst).^2 + imag(delay_prod_inst).^2);
error_mag = mag_double - mag;
% ошибка в процентах
er = error_mag ./ mag_double;
error_in_percent = round(abs(er) .* 100); 

%% Знаменатель
complex_mult = rxWaveform .* conj(rxWaveform);

% avg_channel
complex_mult_avg = complex_mult;

for i = 2 : size(complex_mult)
    if (i < 17)
        complex_mult(i) = complex_mult(i) + complex_mult(i-1);
    else
        complex_mult(i) = complex_mult(i-1) + complex_mult_avg(i) - complex_mult_avg(i-16);
    end
end

complex_mult = complex_mult .* 0.75;

%%
count = 0;
for i = 1:size(mag)
    if (mag(i) > complex_mult(i))
        r(i) = 1;
        count = count + 1;
    else
        r(i) = 0;
    end

    if (count > 100)
        r(i+1) = 0;
        break;
    end
end
plot(r)



%%

% s = size(rxWaveform);
% r = zeros(floor(s(1)/16),1);
% 
%  for i = 1:size(r)
%      aa(i) =  sum(delay_prod_inst((i*16)-15:(i*16)));
%      bb(i) =  sum(e((i*16)-15:(i*16)));
%      c(i) = bb(i)*threshold;
%      if (aa(i) > c(i)) 
%         r(i) = 1;
%      end
%  end

% subplot(2,1,1)
% plot(r)
% zx = xcorr(real(rxWaveform(1:160)), 'coeff');
% subplot(2,1,2)
% plot(zx)


% Генерация тестовых входных файлов для моделирования
generation = 0;

switch generation
    case 1

    waveform_i = floor((real(rxWaveform))*2^14);
    waveform_q = floor((imag(rxWaveform))*2^14);

    writematrix(waveform_i, 'E:\Netlist\ITGLOBAL\wi_fi\correlator_src\tb\waveform_i.txt');
    writematrix(waveform_q, 'E:\Netlist\ITGLOBAL\wi_fi\correlator_src\tb\waveform_q.txt');

    otherwise
end

clear all