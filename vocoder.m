% Austin Gilbert
% B EE 235 Spring 2020
% Final Project: N-Channel Vocoder
% vocoder.m

function y = vocoder(inputfile, outputfile, bandnum, soundon, graphon, noiseon)
% This N-Channel Vocoder modulates a signal into Vocoded Speech
%   inputfile - the signal to modulate with the carrier signal
%   outputfile - file name to use for output
%   bandnum - number of bands/channels
%   soundon - the boolean to play sound
%   graphon - the boolean to create figures
%   noiseon - the boolean to replace carrier signal with white noise

% Description
disp(['Modulating ',inputfile,' to ',outputfile,' through a ',num2str(bandnum),'-Band Vocoder.']);
f = waitbar(0, 'Starting');

% Boundary Frequency Generation for N Bands
xmin = log10(300/165.4+1)/2.1;                  % min log frequency
xmax = log10(6000/165.4+1)/2.1;                 % max log frequency
x = xmin:(xmax-xmin)/bandnum:xmax;
fco = zeros(1,bandnum+1);                       % skeleton of fco vector
for i = 1:bandnum+1
    fco(i) = 165.4*(10^(x(i)*2.1)-1);           % create fco vector
end

% Bandpass Filter
[signal, Fs] = audioread(inputfile);            % read inputfile
t = (0:length(signal)-1)/Fs;                    % time vector creation
y = zeros(length(t), 1);
carrier = rand(size(t));                        % noise creation
waitbar(.33, f, 'Modulating');
for n = 1:length(fco)-1
    [bb, aa] = butter(3,[fco(n), fco(n+1)]/(Fs/2));
    if graphon
        figure(1)                               % GRAPH - Freq Plot
        [H, F] = freqz(bb, aa, 256, Fs);
        plot(F, abs(H), 'DisplayName', ['Band ', num2str(n)]);
        xlabel('Frequency (Hz)');
        ylabel('|H|');
        legend;
        title('Frequency Response (Bandpass)');
        hold on;
    end
    bpsignal = filter(bb, aa, signal);          % bandpass filter
    if graphon
        figure(2)                               % GRAPH - Bandpassed Signal
        plot(t, bpsignal, 'DisplayName', ['Band ', num2str(n)]);
        xlabel('Time (s)');
        ylabel('Amplitude (normalized)');
        legend
        title('Bandpassed Signals');
        hold on;
    end
    
% Envelope Rectification
    rectsignal = abs(bpsignal);                 % full-wave rectification
    [dd, cc] = butter(2, 400/(Fs/2));           % lowpass filter
    if graphon
        waitbar(.67, f, 'Graphing Intermediate Waveforms');
        figure(3)                               % GRAPH - Freq Plot
        [H, F] = freqz(dd, cc, 256, Fs);
        plot(F, abs(H), 'DisplayName', ['Band ', num2str(n)]);
        xlabel('Frequency (Hz)');
        ylabel('|H|');
        legend;
        title('Frequency Response (Lowpass)');
        hold on;
    end
    envelope = filter(dd, cc, rectsignal);
    if graphon
        figure(4)                               % GRAPH - Extracted Envelopes
        subplot(bandnum/2, 2, n)
        plot(t, envelope);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title(['Envelope ', num2str(n)]);
        hold on;
    end
    
% Carrier Modulation / Summation
    if ~noiseon                                 % modulation with carrier
        carrier = sin(2*pi*round((fco(n) + fco(n+1))/2)*t);
        if graphon
            figure(5)                           % GRAPH - Modulated Signals
            subplot(bandnum/2, 2, n)
            plot(t, envelope .* carrier');
            xlabel('Time (s)');
            ylabel('Amplitude');
            title(['Mod Signal ', num2str(n)]);
            hold on;
            
        end
        y = y + (envelope .* carrier');         % summation & modulation
    else
        carriermod = envelope .* carrier';
        if graphon
            figure(5)                           % GRAPH - Modulated Signals
            subplot(bandnum/2, 2, n)
            plot(t, filter(bb, aa, carriermod));
            xlabel('Time (s)');
            ylabel('Amplitude');
            title(['Mod Signal ', num2str(n)]);
            hold on;
        end
        y = y + filter(bb, aa, carriermod);     % summation & bandpass
    end
end
y = y * (max(abs(signal))/max(abs(y)));
audiowrite(outputfile, y, Fs);

if graphon
    waitbar(.84, f, 'Graphing Final Waveforms');
    
    figure(6)                                   % GRAPH - Original/Modulated Waveforms
    subplot(211)
    plot(t, y)
    xlabel('Time (s)');
    ylabel('Amplitude (normalized)');
    title('Modulated Signal');
    subplot(212)
    plot(t, signal)
    xlabel('Time (s)');
    ylabel('Amplitude (normalized)');
    title('Original Signal');
    
    figure(7)                                   % GRAPH - Spectrograms
    subplot(211)
    [s, x, a] = spectrogram(y, hamming(512) ,32, 1024, Fs);
    surf(a, x, 20*log10(abs(s)), 'EdgeColor', 'none');
    colormap(jet); view(0,90);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Modulated Signal Spectrogram');
    subplot(212)
    [s, x, a] = spectrogram(signal, hamming(512) ,32, 1024, Fs);
    surf(a, x, 20*log10(abs(s)), 'EdgeColor', 'none');
    colormap(jet); view(0,90);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Original Signal Spectrogram');
end

waitbar(1, f, 'Finishing')
close(f)

if soundon
    linelength = fprintf('The Original Sound'); % SOUND - Original Sound
    sound(signal, Fs)
    pause(2)
    fprintf(repmat('\b',1,linelength))
    
    linelength = fprintf('The Modulated Sound');% SOUND - Modulated Sound
    sound(y, Fs)
    pause(2)
    fprintf(repmat('\b',1,linelength))
end
end