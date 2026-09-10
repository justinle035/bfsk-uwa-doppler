function baseband_sig = linecoding(bits, N, high_bit)
    % bits: input bit sequence (e.g. [0 1 0 1])
    % N: samples per bit (e.g. 9)
    % high_bit: which bit value should map to the high level (e.g. 1 or 0)

    % Mark logic-1 wherever the bit equals high_bit, logic-0 otherwise
    % This makes the "high" level assignment flexible for either bit value
    logic_bits = (bits == high_bit);

    % Expand into the sample-rate waveform
    baseband_sig = repmat(logic_bits(:), 1, N)';
    baseband_sig = baseband_sig(:)';
end
