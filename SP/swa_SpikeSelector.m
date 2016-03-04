function swa_SpikeSelector(EEG, MONTAGE)
OPTS.window_len = 1; % in seconds
OPTS.min_spike_freq = 7;
OPTS.max_spike_freq = 25;
OPTS.wavelet_freq_spacing = 2; % in Hz
OPTS.wavelet_name = 'mexh';
OPTS.threshold_histogram = 0;

SPIKE_BOUNDS = NaN(length(EEG.csc_event_data), 2);

HANDLES = DefineInterface();

GREY = [0.6 0.6 0.6];
CURR_SPIKE = 1;
CURR_CHAN = 1;
PlotSpike();

    function handles = DefineInterface()
    
    % Dual monitors creates an issue in Linux environments whereby the two
    % screens are seen as one long screen, so always use first one
    sd = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment.getScreenDevices;
    fig_height = 1500; % in pixels
    if length(sd) > 1
        bounds = sd(2).getDefaultConfiguration.getBounds;
        figPos = [bounds.x bounds.y bounds.width fig_height];
    else
        bounds = sd(1).getDefaultConfiguration.getBounds;
        figPos = [bounds.x bounds.y bounds.width fig_height];
    end
        
    handles.selector = figure(...
        'Name', 'Spike Selector', ...
        'color', 'w', ...
        'KeyPressFcn', @cb_KeyPressed, ...
        'Units', 'pixels', ...
        'OuterPosition', figPos);
    handles.rawdata_ax = axes(...
        'Parent', handles.selector, ...
        'Title', 'EEG Time Course', ...
        'xlim', [-OPTS.window_len , OPTS.window_len], ...
        'OuterPosition', [0 .66 1 0.33]);
    handles.wavelet_ax = axes(...
        'Parent', handles.selector, ...
        'Title', 'Mean Wavelet Energy', ...
        'xlim', [-OPTS.window_len, OPTS.window_len], ...
        'OuterPosition', [0 .33 1 .33]);
    handles.hist_ax = axes(...
        'Parent', handles.selector, ...
        'Title', 'Wavelet Peak Energies Distribution', ...
        'OuterPosition', [0 0 1 .33]);
    end

    function cb_KeyPressed(hObject, cbData)
    num_spikes = length(EEG.csc_event_data);
    num_chans = length(MONTAGE.data);
    % movement keys
    switch cbData.Key
        case 'rightarrow'
            % move to the next spike if not at the end
            if  CURR_SPIKE < num_spikes
                CURR_SPIKE = CURR_SPIKE + 1;
                PlotSpike();
            end
        case 'leftarrow'
            if CURR_SPIKE > 1
               CURR_SPIKE = CURR_SPIKE - 1;
                PlotSpike();
            end
        case 'downarrow'
            % move to the next channel in montage if not at the end
            if  CURR_CHAN < num_chans
                HANDLES.rawdata(CURR_CHAN).LineWidth = 1;
                HANDLES.rawdata(CURR_CHAN).Color = GREY;
                CURR_CHAN = CURR_CHAN + 1;
                HANDLES.rawdata(CURR_CHAN).LineWidth = 2;
                HANDLES.rawdata(CURR_CHAN).Color = 'Black';
            end
        case 'uparrow'
            if CURR_CHAN > 1
               HANDLES.rawdata(CURR_CHAN).LineWidth = 1;
               HANDLES.rawdata(CURR_CHAN).Color = GREY;
               CURR_CHAN = CURR_CHAN - 1;
               HANDLES.rawdata(CURR_CHAN).LineWidth = 2;
               HANDLES.rawdata(CURR_CHAN).Color = 'Black';
            end
        case 'x' % save spike bounds to workspace
            assignin('base', 'SPIKE_BOUNDS', SPIKE_BOUNDS);
    end
    end

    function cb_click(hObject, cbData)
        switch HANDLES.selector.SelectionType
        case 'normal' % left click
            SPIKE_BOUNDS(CURR_SPIKE, 1) = cbData.IntersectionPoint(1);
            PlotSpike();
        case 'alt' % right click
            SPIKE_BOUNDS(CURR_SPIKE, 2) = cbData.IntersectionPoint(1);
            PlotSpike();
        end
    end

    function PlotSpike()
    % get spike info
    spiketime = EEG.csc_event_data{CURR_SPIKE, 2};
    spikesample = floor(spiketime .* EEG.srate);

    % get window around spike
    samples_per_window = OPTS.window_len * EEG.srate;
    window = [spikesample - samples_per_window : spikesample + samples_per_window];

    % select the data in this window, as per the montage used to mark spikes
    montage_chans = EEG.data([MONTAGE.data{:, 2}], window);
    montage_refs = EEG.data([MONTAGE.data{:, 3}], window);
    window_data = montage_chans - montage_refs;

    
    % Plot EEG Time Course
    % ^^^^^^^^^^^^^^^^^^^
    % butterfly plot with raw data
    time_range = [-samples_per_window : samples_per_window] / EEG.srate;
    HANDLES.rawdata = plot(HANDLES.rawdata_ax, time_range, window_data, ...
        'color', GREY, 'buttondownfcn', @cb_click);
    
    % highlect selected channel
    HANDLES.rawdata(CURR_CHAN).LineWidth = 2;
    HANDLES.rawdata(CURR_CHAN).Color = 'Black';
    
    % plot vertical bars deliminiting spike start/end
    axes(HANDLES.rawdata_ax);
    spike_start = SPIKE_BOUNDS(CURR_SPIKE, 1);
    spike_end = SPIKE_BOUNDS(CURR_SPIKE, 2);
    if ~isnan(spike_start)
        h = line([spike_start spike_start], HANDLES.rawdata_ax.YLim, ...
            'Color', 'green');
    end
    if ~isnan(spike_end)
        line([spike_end spike_end], HANDLES.rawdata_ax.YLim, ...
            'Color', 'red');
    end
    
    % Plot wavelet
    % ^^^^^^^^^^^^
    axes(HANDLES.wavelet_ax);
    % design wavelets used for spike detection
    wavelet_freqs = OPTS.min_spike_freq : OPTS.wavelet_freq_spacing ...
        : OPTS.max_spike_freq;
    wavelet_scales = centfrq(OPTS.wavelet_name) ./ wavelet_freqs * EEG.srate;
    
    % calculate wavelets of spike channels
    [num_chans, datalen] = size(window_data);
    spike_wavelets = nan(length(wavelet_scales), datalen, num_chans);
    for n_channel = 1 : num_chans
        % find the wavelets corresponding to the pseudo-frequencies
        x = cwt(window_data(n_channel, :), wavelet_scales, ...
            OPTS.wavelet_name);
        % average the wavelets over the scale-range
        x = mean(x, 1);
        % find the largest peak
        [peak_amp(n_channel), peak_ind(n_channel)] = max(abs(x));
    end
    plot(HANDLES.wavelet_ax, time_range, x);
    
    % Plot Histogram 
    % ^^^^^^^^^^^^^^
    axes(HANDLES.hist_ax);
    if OPTS.threshold_histogram
        hist(peak_ind(peak_amp > prctile(peak_amp, 80)));
    else
        hist(peak_ind);
    end
    
    % potential range
    % potential_range = peak_ind(peak_amp > prctile(peak_amp, 80)) - OPTS.window_len;
    
    end
end


