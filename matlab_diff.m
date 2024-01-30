% Check your txt file against audio file in Matlab or Python

% Original Audio Wav File
y = audioread('/Users/nicolettecash/ase_projects/ase-2024/sweep.wav');

% Rust Text File
file = fopen('text.txt','r');
formatting = '%f';
z = fscanf(file,formatting);

col1 = z(1:132300);
col2 = z(132301:264600);

rust_audio = [col1 col2];

% Calculate difference between rust and matlab
diff = y - rust_audio;
