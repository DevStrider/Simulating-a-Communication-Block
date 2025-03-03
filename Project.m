a = 1; 
b = 8; 
t = 0:0.3:30;
x_t = a * sin(0.5*pi*t*b);
    
% Sampling Task 1
Ns = 2; 
ts = t(1 : Ns : end);
xs = x_t(1 : Ns : end);

%Quantization Task 2 
function xq = uniform_quantizer(xs , L)
vmax = max(xs);
vmin = min(xs);
delta = (vmax - vmin)/L;
quantized_levels = vmin+(delta/2) : delta : vmax-(delta/2);
for i = 1 : length(xs)
    abs_error = abs(xs(i) - quantized_levels);
    l = find(abs_error == min(abs_error));
    xq(i) = quantized_levels(l(1));
end
end

%Plot Mean Absolute Quantization Error vs Number of Levels Task 3
L = [2 4 8];
mean_abs_error = zeros(size(L));
for i = 1 : length(L)
    xq = uniform_quantizer(xs , L(i));
    mean_abs_error(i) = mean(abs(xs - xq));
end
figure; 
plot(L , mean_abs_error, '-o');
xlabel('Number of Quantization Levels (L)');
ylabel('Mean Absolute Quantization Error');
title('MAQ Error vs Number of Quantization Levels');
grid on;

% Plot Practical and Theoretical Varience vs Number of levels Task 4
pratical_varience = zeros(size(L));
theoretical_variance = zeros(size(L));
for i = 1:length(L)
    xq = uniform_quantizer(xs , L(i));
    pratical_varience(i) = var(xs - xq);
    vmax = max(xs);
    vmin = min(xs);
    delta = (vmax - vmin)/L(i);
    theoretical_variance(i) = (delta^2)/12;
end
figure;
plot(L, pratical_varience, '-o', 'DisplayName', 'Practical Variance');
hold on;
plot(L, theoretical_variance, '-s', 'DisplayName', 'Theoretical Variance');
xlabel('Number of Quantization Levels (L)');
ylabel('Variance of Quantization Error');
title('Quantization Error Variance vs Number of Levels');
legend('Location', 'best');
grid on;
hold off;

% Plot Practical and Theoretical SQNR vs Number of levels Task 5
practical_SQNR = zeros(size(L));
theoretical_SQNR = zeros(size(L));
for i = 1:length(L)
    practical_SQNR(i) = (a^2) / pratical_varience(i);
    theoretical_SQNR(i) = 3 * L(i)^2;
end
%Convert SQNR into dB
% practical_SQNR_dB = 10 * log10(practical_SQNR); 
% theoretical_SQNR_dB = 10 * log10(theoretical_SQNR);
figure;
plot(L, practical_SQNR, '-o', 'DisplayName', 'Practical SQNR');
hold on;
plot(L, theoretical_SQNR, '-s', 'DisplayName', 'Theoretical SQNR');
xlabel('Number of Quantization Levels (L)');
ylabel('SQNR');
title('SQNR vs Number of Quantization Levels');
legend('Location', 'best');
grid on;
hold off;

L = 8;
xq = uniform_quantizer(xs , L);
delta = (vmax - vmin)/L;

% Encode quantized levels Task 6
uniqueQuantizedLevels = vmin+(delta/2) : delta : vmax-(delta/2);%Unique quantized values in ascending order
encoded = zeros(1,length(xq));
for i = 1:length(xq)
    for j = 1:length(uniqueQuantizedLevels)
        if xq(i) == uniqueQuantizedLevels(j)
            encoded(i) = j-1;
        end
    end
end

% Huffman encoder Task 7
unique_values = unique(encoded);
[counts, edges] = histcounts(encoded, [unique_values - 0.5, unique_values(end) + 0.5]);
probabilities = counts / length(encoded);
[dictionary , x] = huffmandict(unique(encoded),probabilities);
sourceEncoded = huffmanenco(encoded,dictionary);

% Noiseless channel Task 8
noise = 0;
sourceEncoded = sourceEncoded + noise;

% Huffman decoder Task 10
sourceDecoded = huffmandeco(sourceEncoded,dictionary);

% Decoder Task 9
uniqueDecodedLevels = zeros(1,L);
uniqueDecodedLevels(1) = vmin + (delta/2);
uniqueDecodedLevels(L) = vmax - (delta/2);
for i=2:(L-1)
    uniqueDecodedLevels(i) = uniqueDecodedLevels(i-1) + delta;
end
decodedSequence = zeros(1,length(sourceDecoded));
for i = 1:length(sourceDecoded)
    decodedSequence(i) = uniqueDecodedLevels(sourceDecoded(i) + 1);
end

% Plot Input and Output Signals Task 11
figure;
plot(t, x_t, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Input Signal (x_t)');
hold on;
plot(ts, decodedSequence, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Output Signal (Decoded)');
xlabel('Time (s)');
ylabel('Amplitude');
title('Input vs Output Signal');
legend('Location', 'best');
grid on;
hold off;

% Efficiency and CR Task 13 and 14
I = -log2(probabilities);
H = sum(probabilities .* I);
L = sum(probabilities .* ceil(I));
eff = H / L;
CR = log2(length(unique(encoded))) / L;
disp(['Compression Efficiency: ', num2str(eff)]);
disp(['Compression Rate (CR): ', num2str(CR)]);