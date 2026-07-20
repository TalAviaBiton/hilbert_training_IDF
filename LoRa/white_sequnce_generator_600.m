reg = uint8(255); % 0xFF
seq = zeros(1, 600, 'uint8');
for i = 1:600
    seq(i) = reg;
    % חישוב הביט הבא לפי הפולינום x^8+x^6+x^5+x^4+1
    new_bit = xor(bitget(reg, 8), xor(bitget(reg, 6), xor(bitget(reg, 5), bitget(reg, 4))));
    reg = bitxor(int32(bitshift(reg, 1)), int32(new_bit));
end
% הדפסת הרצף למסך כדי להעתיק לקוד
fprintf('uint8([%s])\n', num2str(seq, '%d, '));