clc;
clear;
close all;

%% Function Definitions

% Function to calculate DNL and INL
function [dnl_result, inl_result, avg_dnl_result, avg_inl_result] = compute_dnl_inl(adc_codes, code_hits, full_scale, amplitude, bits)
    % Remove NaN values from the data
    adc_codes = adc_codes(~isnan(adc_codes));
    code_hits = code_hits(~isnan(code_hits));

    % Total number of hits across all codes
    total_samples = sum(code_hits);

    % Compute the theoretical probability density function for each code
    prob_density = (1 / pi) * ( ...
        asin((full_scale * (adc_codes - 2^(bits-1)) / (amplitude * 2^bits))) - ...
        asin((full_scale * (adc_codes - 1 - 2^(bits-1)) / (amplitude * 2^bits))) ...
    );

    % Differential Nonlinearity (DNL) calculation
    dnl_result = (code_hits ./ (total_samples * prob_density)) - 1;

    % Ensure DNL values are real
    dnl_result = real(dnl_result);

    % Integral Nonlinearity (INL) calculation
    inl_result = cumsum(dnl_result);

    % Ensure INL values are real
    inl_result = real(inl_result);

    % Compute average DNL and INL values
    avg_dnl_result = mean(dnl_result);
    avg_inl_result = mean(inl_result);
end


% Function to plot histogram, DNL, and INL separately
function plot_histo_INL_DNL(frequency, adc_codes, code_hits, dnl_result, inl_result, plot_title)
    % Plot Histogram
    figure('Position', [200, 200, 800, 600]);
    bar(adc_codes, code_hits, 'FaceColor', [0.4, 0.8, 0.4], 'EdgeColor', 'k'); % Green bars
    xlabel('ADC Code', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Hits', 'FontSize', 12, 'FontWeight', 'bold');
    title([plot_title, ' - Histogram'], 'FontSize', 14, 'FontWeight', 'bold');
    grid on;

    % Plot DNL
    figure('Position', [300, 300, 800, 600]);
    plot(adc_codes, dnl_result, 'b-', 'LineWidth', 1.5); % Dashed blue line
    xlabel('ADC Code', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('DNL (LSBs)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Differential Nonlinearity (DNL) - Frequency: %.2f Hz', frequency), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    xlim([min(adc_codes), max(adc_codes)]);

    % Plot INL
    figure('Position', [400, 400, 800, 600]);
    plot(adc_codes, inl_result, 'b-', 'LineWidth', 1.5); % Solid blue line
    xlabel('ADC Code', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('INL (LSBs)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Integral Nonlinearity (INL) - Frequency: %.2f Hz', frequency), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    xlim([min(adc_codes), max(adc_codes)]);
end


%% Main Calculation and Plotting

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
        new_codes(1, k) = adc_decoded(index);
    end

    % Create histogram with step size of 1
    codes = min(new_codes):1:max(new_codes);
    hits = histcounts(new_codes, [codes, max(codes) + 1]);

    FSR = 2.5;               % Input voltage
    Amp = FSR / 2;           % Sine wave amplitude
    resolution = 9;          % Adc resolution 

    % Calculate DNL and INL
    [dnl, inl, avg_dnl, avg_inl] = compute_dnl_inl(codes, hits, FSR, Amp, resolution);

    % Display average DNL and INL
    fprintf('Frequency: %.2f Hz - Average DNL: %.2f, Average INL: %.2f\n', current_frequency, avg_dnl, avg_inl);

    % Plot combined graph
    plot_histo_INL_DNL(current_frequency, codes, hits, dnl, inl, sprintf('Frequency at %.2f Hz', current_frequency));
end
