% ===== Parameters (EDIT THESE) =====
fs        = 44100;   % Sampling rate (Hz)
duration  =20;       % Duration (seconds)

f_carrier = 350;     % Carrier frequency (Hz)
f_mod     = 40;      % Modulation frequency (Hz)

amp_min   = 0.2;     % Minimum amplitude
amp_max   = 0.9;     % Maximum amplitude

volume    = 0.4;     % 🔊 Master volume (0 → 1)

% ===== Time vector =====
t = 0:1/fs:duration;

% ===== Modulator (scaled sine from amp_min to amp_max) =====
modulator = (sin(2*pi*f_mod*t) + 1)/2;          % 0 → 1
modulator = amp_min + (amp_max - amp_min)*modulator;

% ===== Carrier signal =====
carrier = sin(2*pi*f_carrier*t);

% ===== Amplitude-modulated signal =====
signal = modulator .* carrier;

% ===== Apply volume =====
signal = volume * signal;

% ===== Play sound =====
sound(signal, fs);

%% ===== Optional: visualize =====
figure;
subplot(3,1,1); plot(t, modulator); title('Modulator (40 Hz)');
subplot(3,1,2); plot(t, carrier); title('Carrier (350 Hz)');
subplot(3,1,3); plot(t, signal); title('AM Signal');
xlabel('Time (s)');