function txt = bits2text_pad(bits)
%BITS2TEXT_PAD  Convert a 1xN bit vector back to an ASCII string
% Automatically zero-pads the tail if the bit count is not a multiple of 8

    bits = bits(:).';   % ensure row vector

    r = mod(length(bits),8);
    if r ~= 0
        % zero-pad the tail
        bits = [bits zeros(1, 8-r)];
    end

    % reshape into bytes
    b = reshape(bits, 8, []).';

    % bits -> number
    ascii = bi2de(b, 'left-msb');

    % number -> character
    txt = char(ascii).';
end
