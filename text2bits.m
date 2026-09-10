function bits = text2bits(txt)
%TEXT2BITS  Convert a character string to a bit vector (1xN)
%   txt  : char array or string
%   bits : 1xN vector of 0/1 bits (MSB first)

    if isstring(txt)
        txt = char(txt);
    end

    % Convert characters to ASCII codes
    ascii = uint8(txt);

    % Each character -> 8 bits, MSB first
    b = de2bi(ascii, 8, 'left-msb');

    % Flatten into a single row vector
    bits = reshape(b.', 1, []);
end
