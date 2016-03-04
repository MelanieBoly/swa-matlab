adjustment_window = 500; % in samples
time_range = [- adjustment_window : adjustment_window] ...
    / EEG.srate;

% get 10 - 20 montage
montage = load('10 - 20.emo', '-mat');

% wavelet 
wavelet_name = 'mexh';
freq_range = 8 : 2 : 16;
spike_scales = centfrq(wavelet_name)./ freq_range * EEG.srate;


% Load cleaned, annotated data.
spikefile = uigetfullfile('.set');
EEG = pop_loadset(spikefile);

% Automatically detect spike peaks based on manual annotations.
% FIXME: Assuming all events are spikes
for n = 1 : length(EEG.csc_event_data)
    spike_time = EEG.csc_event_data{n, 2}; % in seconds
    spike_sample = floor(spike_time * EEG.srate);
    adjustment_window_range = spike_sample - adjustment_window ...
        : spike_sample + adjustment_window;
    
    % calculate the spike data using the montage channels
    spike_channel = EEG.data([montage.data{:, 2}], adjustment_window_range);
    spike_reference = EEG.data([montage.data{:, 3}], adjustment_window_range);
    spike_data = spike_channel - spike_reference;
          
    if flag_plot
       figure('color', 'w');
       axes('nextPlot', 'add', ...
           'xlim', [time_range(1) , time_range(end)]);
       plot(time_range, ...
           spike_data, ...
           'color', [0.6, 0.6, 0.6]);
    end
    
    % calculate wavelets of spike channels
    spike_wavelets = nan(length(spike_scales), size(spike_data, 2), size(spike_data, 1));
    for n_channel = 1 : size(spike_data, 1)
        % find the wavelets corresponding to the pseudo-frequencies
        x = cwt(spike_data(n_channel, :), spike_scales, wavelet_name);
        % average the wavelets over the scale-range
        x = mean(x, 1);
        % find the largest peak
        [peak_amp(n_channel), peak_ind(n_channel)] = max(abs(x));
    end
    
     % potential range
     potential_range = peak_ind(peak_amp > prctile(peak_amp, 80)) - adjustment_window;
    
end
   


