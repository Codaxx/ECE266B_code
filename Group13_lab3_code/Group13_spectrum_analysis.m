clc;
clear;
close all;

% Load and process frequency array
frequency_list = readtable('sweeping_frequencies.csv');
frequency_list = floor(table2array(frequency_list));

% Sort the frequencies in ascending order and remove duplicates
frequency_list = unique(sort(frequency_list));

% Prepare storage for results
results = zeros(length(frequency_list), 4); % Columns: [Frequency, SNR, SNDR, SFDR]

% Loop through each frequency
for idx = 1:length(frequency_list)
    current_frequency = frequency_list(idx);
    results(idx, 1) = current_frequency;

    % Generate filename and read data
    filepath = ''; % Update if needed
    filename = strcat(filepath, 'INPUT_SIN_', num2str(current_frequency), '.0Hz_.csv');
    raw_data = table2array(readtable(filename));

    % Process data
    binary_data = dec2bin(raw_data);
    end_of_conversion = str2num(binary_data(:, 1));
    adc_codes = str2num(binary_data(:, 2:10));
    adc_decoded = bin2dec(num2str(adc_codes));
    new_codes = size(length(adc_decoded), 1);

    % Find rising edges and extract new codes
    rising_edges = find(diff(end_of_conversion) == 1) + 1;
    for k = 1:length(rising_edges)
        index = rising_edges(k);
        new_codes(k, 1) = adc_decoded(index);
    end

    % Analyze data and calculate SNR, SNDR, SFDR
    sampling_rate = 1e6;
    num_segments = 1;
    bandwidth = sampling_rate / 2;
    signal_data = new_codes;
    periodogram_length = floor (length(signal_data) / num_segments);
    periodogram_length = round(periodogram_length);
    % plot_settings = [0, 0, 1, 0, 1]; % [plotYN, plotAll, plotHold, plotLin, datNorm]
    
    [sinusoid_power, data_minus_sinusoid_in_BW_power, SNDR, ENOB,HD2,HD3, SNR, SFDR] = ...
        plot_periodogram_SFDR(signal_data, periodogram_length, num_segments, current_frequency, ...
        sampling_rate, bandwidth, 0, 0, 1, 0, 1);
    
    fprintf('Frequency = %.3f, SNR = %.3f dB, SNDR = %.3f dB, SFDR = %.3f dB, HD2 = %.3f, HD3 = %.3f\n', ...
        current_frequency, SNR, SNDR, SFDR, HD2, HD3);
    
    results(idx, 2:4) = [SNR, SNDR, SFDR];

    % Reconstruct and plot waveform
    adc_resolution = 9;
    v_min = 0; v_max = 1;
    mapped_voltages = ((signal_data / (2^adc_resolution - 1)) * (v_max - v_min)) + v_min;
    time_axis = (0:length(signal_data) - 1) / sampling_rate * 1000; % Time in ms
    
    figure;
    plot(time_axis, mapped_voltages, 'b');
    title('Reconstructed Time-Domain Waveform');
    xlabel('Time (ms)');
    ylabel('Amplitude (V)');
    grid on;
end

% Plot SNR, SNDR, SFDR results
figure;
semilogx(results(:, 1), results(:, 2), 'b-', 'LineWidth', 1.5, "DisplayName", "SNR");
hold on;
semilogx(results(:, 1), results(:, 3), 'm-', 'LineWidth', 1.5, "DisplayName", "SNDR");
semilogx(results(:, 1), results(:, 4), 'r-', 'LineWidth', 1.5, "DisplayName", "SFDR");
hold off;
xlabel("Frequency (Hz)");
ylabel("SNR / SNDR / SFDR (dB)");
title("SNR, SNDR, and SFDR vs. Frequency");
legend('Location', 'best');
grid on;

