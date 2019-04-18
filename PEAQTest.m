function odg = PEAQTest(filename, fdir)

disp(filename);

orig = [fdir filename '.wav'];
ref = [fdir filename '_48k.wav'];
test = [fdir 'decoded/' filename '_decoded.wav'];

% resample the original file to 48k sr
if exist(ref, 'file') ~= 2
    [y,fs] = audioread(orig);
    y = resample(y, 48000, fs);
    audiowrite(ref, y, 48000);
end

decodedFile = full_codec_separate(ref, 64000, test);

[odg, movb] = PQevalAudio_fn(ref, test);