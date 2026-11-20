// ============================================================================
// EXACT bounded-posit design. Arithmetic modules taken from the
// verified bmac_exact.v (includes product saturation). bmac, bsingalPE and
// b3_3systolic are unchanged from the original design except the multiplier
// instance (posit_mult_exact).
// Read ONLY this file in a synthesis run (module names match the other variant).
// ============================================================================
// ============================================================================
// Bounded-posit decoder
// ============================================================================
module bposit_decode #(
    parameter integer N  = 12,
    parameter integer es = 2,
    parameter integer rS = 4,
    parameter integer Bs =
        (N <= 8)  ? 3 :
        (N <= 16) ? 4 :
        (N <= 32) ? 5 :
        (N <= 64) ? 6 : 7
)(
    input  wire [N-1:0]     in,
    output wire             sign,
    output wire             is_zero,
    output wire             is_nar,
    output wire             rc,
    output wire [Bs:0]      regime,
    output wire [es-1:0]    exp,
    output wire [N-es-2:0]  frac,
    output wire             frac_zero
);

    localparam integer MUXW  = N - 3;
    localparam integer MW    = N - es - 1;
    localparam integer NCAND = rS - 1;

    assign sign = in[N-1];

    wire [N-2:0] body = in[N-2:0];
    assign rc = in[N-2];

    wire body_nonzero = |body;
    assign is_nar  =  sign & ~body_nonzero;
    assign is_zero = ~sign & ~body_nonzero;

    wire [rS-2:0] xbit;
    wire [rS-1:0] sel;

    genvar gi;
    generate
        for (gi = 1; gi < rS; gi = gi + 1) begin : GEN_XBIT
            assign xbit[gi-1] = in[N-2-gi] ^ rc;
        end
    endgenerate

    assign sel[0] = xbit[0];
    generate
        for (gi = 2; gi < rS; gi = gi + 1) begin : GEN_SEL
            assign sel[gi-1] = xbit[gi-1] & ~(|xbit[gi-2:0]);
        end
    endgenerate
    assign sel[rS-1] = ~(|xbit);

    reg [Bs:0] run_len;
    integer k;
    always @(*) begin
        run_len = {(Bs+1){1'b0}};
        for (k = 1; k <= rS; k = k + 1)
            if (sel[k-1])
                run_len = k;
    end

    assign regime = rc
        ? (run_len - {{Bs{1'b0}}, 1'b1})
        : run_len;

    wire [NCAND*MUXW-1:0] cand_flat;

    generate
        for (gi = 2; gi <= rS; gi = gi + 1) begin : GEN_CAND
            if (gi == 2) begin : GEN_WIDEST_CAND
                assign cand_flat[(gi-2)*MUXW +: MUXW] =
                    body[N-2-gi -: (N-1-gi)];
            end else begin : GEN_PADDED_CAND
                assign cand_flat[(gi-2)*MUXW +: MUXW] = {
                    body[N-2-gi -: (N-1-gi)],
                    {(gi-2){1'b0}}
                };
            end
        end
    endgenerate

    reg [MUXW-1:0] exp_sig;
    always @(*) begin
        exp_sig = cand_flat[0 +: MUXW];
        for (k = 1; k < rS; k = k + 1)
            if (sel[k-1])
                exp_sig = cand_flat[(k-1)*MUXW +: MUXW];

        if (sel[rS-1])
            exp_sig = cand_flat[(rS-2)*MUXW +: MUXW];
    end

    wire [es-1:0] raw_exp;
    wire [MW-1:0] raw_frac;

    assign raw_exp  = exp_sig[MUXW-1 -: es];
    assign raw_frac = {
        exp_sig[MUXW-es-1:0],
        {(MW-(MUXW-es)){1'b0}}
    };

    assign exp       = raw_exp;
    assign frac      = raw_frac;
    assign frac_zero = ~|raw_frac;

endmodule

// ============================================================================
// Bounded scale splitter
// ============================================================================
module bposit_reg_exp_op #(
    parameter integer es = 2,
    parameter integer Bs = 4,
    parameter integer rS = 4
)(
    input  wire signed [Bs+es+1:0] scale_in,
    output reg                      rc,
    output reg  [es-1:0]           exp_o,
    output reg  [Bs:0]             reg_o
);

    localparam [Bs:0] RS_VALUE    = rS;
    localparam [Bs:0] POS_REG_MAX = rS - 1;

    wire signed [Bs+1:0] k = scale_in >>> es;
    wire [es-1:0] e = scale_in[es-1:0];
    wire signed [Bs+1:0] positive_limit = {1'b0, POS_REG_MAX};
    wire [Bs+1:0] negative_limit = {1'b0, RS_VALUE};

    reg [Bs+1:0] k_magnitude;

    always @(*) begin
        k_magnitude = (~k) + {{(Bs+1){1'b0}}, 1'b1};

        if (k >= 0) begin
            rc = 1'b1;
            if (k > positive_limit) begin
                reg_o = POS_REG_MAX;
                exp_o = {es{1'b1}};
            end else begin
                reg_o = k[Bs:0];
                exp_o = e;
            end
        end else begin
            rc = 1'b0;
            if (k_magnitude > negative_limit) begin
                reg_o = RS_VALUE;
                exp_o = {es{1'b0}};
            end else begin
                reg_o = k_magnitude[Bs:0];
                exp_o = e;
            end
        end
    end

endmodule

// ============================================================================
// Bounded-posit encoder with RNE
// ============================================================================
module bposit_encode #(
    parameter integer N  = 12,
    parameter integer es = 2,
    parameter integer rS = 4,
    parameter integer Bs =
        (N <= 8)  ? 3 :
        (N <= 16) ? 4 :
        (N <= 32) ? 5 :
        (N <= 64) ? 6 : 7
)(
    input  wire                    sign,
    input  wire signed [Bs+es+1:0] scale,
    input  wire [N-es:0]           mant_normalized,
    input  wire                    sticky_in,
    input  wire                    special_nar,
    input  wire                    special_zero,
    output wire [N-1:0]            out
);

    localparam integer BODYW   = N - 1;
    localparam integer MFRAC_W = N - es;
    localparam integer NCAND   = rS - 1;
    localparam [Bs:0] RS_VALUE = rS;

    wire rc_o;
    wire [es-1:0] exp_o;
    wire [Bs:0] reg_o;

    bposit_reg_exp_op #(
        .es(es), .Bs(Bs), .rS(rS)
    ) split_scale (
        .scale_in(scale),
        .rc(rc_o),
        .exp_o(exp_o),
        .reg_o(reg_o)
    );

    wire [Bs:0] regime_run = rc_o
        ? (reg_o + {{Bs{1'b0}}, 1'b1})
        : reg_o;

    wire regime_capped = (regime_run >= RS_VALUE);
    wire [Bs:0] regime_width = regime_capped
        ? RS_VALUE
        : (regime_run + {{Bs{1'b0}}, 1'b1});

    wire [MFRAC_W-1:0] mfrac = mant_normalized[MFRAC_W-1:0];

    wire [NCAND*BODYW-1:0] body_candidates;
    wire [NCAND-1:0] guard_candidates;
    wire [NCAND-1:0] round_candidates;
    wire [NCAND-1:0] sticky_candidates;

    genvar wi;
    generate
        for (wi = 2; wi <= rS; wi = wi + 1) begin : GEN_PACK
            wire [wi-1:0] regime_field;
            wire [BODYW-1:0] regime_full;
            wire [BODYW-1:0] exponent_full;
            wire [N-2-wi-es:0] fraction_bits;
            wire [BODYW-1:0] fraction_full;

            if (wi == rS) begin : GEN_MAX_WIDTH_REGIME
                assign regime_field = rc_o
                    ? (regime_capped ? {wi{1'b1}} : {{(wi-1){1'b1}}, 1'b0})
                    : (regime_capped ? {wi{1'b0}} : {{(wi-1){1'b0}}, 1'b1});
            end else begin : GEN_NORMAL_WIDTH_REGIME
                assign regime_field = rc_o
                    ? {{(wi-1){1'b1}}, 1'b0}
                    : {{(wi-1){1'b0}}, 1'b1};
            end

            assign regime_full =
                {{(BODYW-wi){1'b0}}, regime_field} << (BODYW-wi);

            assign exponent_full =
                {{(BODYW-es){1'b0}}, exp_o} << (N-1-wi-es);

            assign fraction_bits =
                mfrac[MFRAC_W-1 -: (N-1-wi-es)];

            assign fraction_full =
                {{(wi+es){1'b0}}, fraction_bits};

            assign body_candidates[(wi-2)*BODYW +: BODYW] =
                regime_full | exponent_full | fraction_full;

            assign guard_candidates[wi-2]  = mfrac[wi];
            assign round_candidates[wi-2]  = mfrac[wi-1];
            assign sticky_candidates[wi-2] =
                (|mfrac[wi-2:0]) | sticky_in;
        end
    endgenerate

    reg [BODYW-1:0] body_pre;
    reg guard_bit;
    reg round_bit;
    reg sticky_bit;
    integer selected_width;

    always @(*) begin
        body_pre   = body_candidates[0 +: BODYW];
        guard_bit  = guard_candidates[0];
        round_bit  = round_candidates[0];
        sticky_bit = sticky_candidates[0];

        for (selected_width = 2;
             selected_width <= rS;
             selected_width = selected_width + 1) begin
            if (regime_width == selected_width) begin
                body_pre = body_candidates[
                    (selected_width-2)*BODYW +: BODYW
                ];
                guard_bit  = guard_candidates[selected_width-2];
                round_bit  = round_candidates[selected_width-2];
                sticky_bit = sticky_candidates[selected_width-2];
            end
        end
    end

    wire last_kept_bit = body_pre[0];
    wire round_increment =
          (guard_bit & (round_bit | sticky_bit))
        | (last_kept_bit & guard_bit & ~(round_bit | sticky_bit));

    wire [N-1:0] rounded_sum =
        {1'b0, body_pre} + {{(N-1){1'b0}}, round_increment};

    wire [BODYW-1:0] rounded_body = rounded_sum[N-1]
        ? {BODYW{1'b1}}
        : rounded_sum[BODYW-1:0];

    wire [N-1:0] positive_word = {1'b0, rounded_body};
    wire [N-1:0] negative_word =
        ~positive_word + {{(N-1){1'b0}}, 1'b1};

    assign out = (special_nar | special_zero)
        ? {special_nar, {(N-1){1'b0}}}
        : (sign ? negative_word : positive_word);

endmodule

// ============================================================================
// EXACT bounded-posit multiplier (baseline for the approximate design)
//
// Drop-in replacement for posit_mult_approx: identical ports, parameters,
// decoder (bposit_decode), encoder (bposit_encode), sign, zero/NaR and scale
// logic. The ONLY difference is the significand path:
//
//   approx : M = 1 + f1 + f2            (fraction cross-product removed)
//   exact  : M = (1 + f1) * (1 + f2)    (full significand product)
//
// Because M_exact lies in [1, 4), normalization is a 1-bit right shift, and
// the bits shifted below the encoder's mantissa width are folded into the
// sticky bit so that RNE rounding in bposit_encode is correct.
//
// Requires the modules bposit_decode and bposit_encode from the existing
// design file (compile both files together).
// ============================================================================
module posit_mult_exact #(
    parameter integer N  = 12,
    parameter integer es = 2,
    parameter integer rS = 4,
    parameter integer Bs =
        (N <= 8)  ? 3 :
        (N <= 16) ? 4 :
        (N <= 32) ? 5 :
        (N <= 64) ? 6 : 7
)(
    input  wire [N-1:0] in1,
    input  wire [N-1:0] in2,
    input  wire         start,
    output wire [N-1:0] out,
    output wire         nar,
    output wire         zero,
    output wire         done
);

    // F  : decoded fraction width
    // PW : width of the full significand product (value in [1,4), 2F frac bits)
    // MW : width of the mantissa expected by bposit_encode (1 + (N-es) bits)
    localparam integer F  = N - es - 1;
    localparam integer PW = 2 * F + 2;
    localparam integer MW = N - es + 1;   // = F + 2

`ifndef SYNTHESIS
    initial begin
        if (N < 4)       $display("ERROR: N must be at least 4");
        if (es < 1)      $display("ERROR: es must be at least 1");
        if (rS < 2)      $display("ERROR: rS must be at least 2");
        if (rS > N-2)    $display("ERROR: rS must be <= N-2");
        if (N < rS+es+2) $display("ERROR: require N >= rS + es + 2");
        if (F < 2)       $display("ERROR: require N - es - 1 >= 2");
        if ((1 << Bs) < N)
            $display("ERROR: Bs is too small for N");
    end
`endif

    // ------------------------------------------------------------------
    // Special values, sign, magnitude (identical to approximate design)
    // ------------------------------------------------------------------
    wire s1 = in1[N-1];
    wire s2 = in2[N-1];

    wire body_nonzero1 = |in1[N-2:0];
    wire body_nonzero2 = |in2[N-2:0];

    wire nar1  =  s1 & ~body_nonzero1;
    wire nar2  =  s2 & ~body_nonzero2;
    wire zero1 = ~s1 & ~body_nonzero1;
    wire zero2 = ~s2 & ~body_nonzero2;

    assign nar  = nar1 | nar2;
    assign zero = (zero1 | zero2) & ~nar;
    assign done = start;

    wire [N-1:0] xin1 = s1 ? (~in1 + {{(N-1){1'b0}}, 1'b1}) : in1;
    wire [N-1:0] xin2 = s2 ? (~in2 + {{(N-1){1'b0}}, 1'b1}) : in2;

    wire rc1, rc2;
    wire [Bs:0]     regime1, regime2;
    wire [es-1:0]   e1, e2;
    wire [F-1:0]    frac1, frac2;
    wire fz1, fz2, zz1, zz2;

    bposit_decode #(.N(N), .es(es), .rS(rS), .Bs(Bs)) decode1 (
        .in(xin1), .sign(), .is_zero(zz1), .is_nar(),
        .rc(rc1), .regime(regime1), .exp(e1),
        .frac(frac1), .frac_zero(fz1)
    );

    bposit_decode #(.N(N), .es(es), .rS(rS), .Bs(Bs)) decode2 (
        .in(xin2), .sign(), .is_zero(zz2), .is_nar(),
        .rc(rc2), .regime(regime2), .exp(e2),
        .frac(frac2), .frac_zero(fz2)
    );

    wire mult_sign = s1 ^ s2;

    // ------------------------------------------------------------------
    // EXACT significand product: (1.f1) x (1.f2), value in [1, 4)
    // ------------------------------------------------------------------
    wire [F:0]    sig1 = {1'b1, frac1};
    wire [F:0]    sig2 = {1'b1, frac2};
    wire [PW-1:0] prod = sig1 * sig2;

    wire mult_m_ovf = prod[PW-1];            // product >= 2

    // Take the leading 1 plus (MW-1) fraction bits; remaining bits -> sticky.
    wire [MW-1:0] mant_hi = prod[PW-1 -: MW];   // used when product >= 2
    wire [MW-1:0] mant_lo = prod[PW-2 -: MW];   // used when product <  2

    wire sticky_hi = |prod[PW-1-MW:0];
    wire sticky_lo = |prod[PW-2-MW:0];

    wire [MW-1:0] mant_normalized = mult_m_ovf ? mant_hi : mant_lo;
    wire          prod_sticky     = mult_m_ovf ? sticky_hi : sticky_lo;

    // ------------------------------------------------------------------
    // Scale combination (identical to approximate design)
    // ------------------------------------------------------------------
    wire signed [Bs+1:0] regime1_ext = {1'b0, regime1};
    wire signed [Bs+1:0] regime2_ext = {1'b0, regime2};

    wire signed [Bs+1:0] r1_signed = rc1 ? regime1_ext : -regime1_ext;
    wire signed [Bs+1:0] r2_signed = rc2 ? regime2_ext : -regime2_ext;

    wire signed [Bs+es+1:0] scale1 = {r1_signed, e1};
    wire signed [Bs+es+1:0] scale2 = {r2_signed, e2};

    wire signed [Bs+es+1:0] scale_result =
        scale1 + scale2 + {{(Bs+es+1){1'b0}}, mult_m_ovf};

    wire [N-1:0] enc_raw;   // encoder output before saturation

    bposit_encode #(.N(N), .es(es), .rS(rS), .Bs(Bs)) encode_result (
        .sign(mult_sign),
        .scale(scale_result),
        .mant_normalized(mant_normalized),
        .sticky_in(prod_sticky),
        .special_nar(nar),
        .special_zero(zero),
        .out(enc_raw)
    );

    // ------------------------------------------------------------------
    // Saturation (same policy as bposit_add): out-of-range products clamp
    // to +/-maxpos or +/-minpos and never flush to zero.
    // ------------------------------------------------------------------
    localparam signed [Bs+es+1:0] FORMAT_MAX_SCALE = (rS * (1 << es)) - 1;
    localparam signed [Bs+es+1:0] FORMAT_MIN_SCALE = -(rS * (1 << es));
    localparam [N-1:0] POS_MAX = {1'b0, {(N-1){1'b1}}};
    localparam [N-1:0] POS_MIN = {{(N-1){1'b0}}, 1'b1};

    wire mult_overflow  = scale_result > FORMAT_MAX_SCALE;
    wire mult_underflow = scale_result < FORMAT_MIN_SCALE;
    wire enc_collapsed  = (enc_raw == {N{1'b0}});
    wire [N-1:0] sat_max = mult_sign ? (~POS_MAX + 1'b1) : POS_MAX;
    wire [N-1:0] sat_min = mult_sign ? (~POS_MIN + 1'b1) : POS_MIN;

    assign out = (nar | zero)    ? enc_raw :
                 mult_overflow   ? sat_max :
                 mult_underflow  ? sat_min :
                 enc_collapsed   ? sat_min :
                                   enc_raw;


endmodule

// ============================================================================
// Standalone bounded-posit adder
// Corrected SIGW = N-es+1
// ============================================================================
module bposit_add #(
    parameter integer N  = 12,
    parameter integer es = 2,
    parameter integer rS = 4,
    parameter integer Bs =
        (N <= 8)  ? 3 :
        (N <= 16) ? 4 :
        (N <= 32) ? 5 :
        (N <= 64) ? 6 : 7
)(
    input  wire [N-1:0] in_a,
    input  wire [N-1:0] in_b,
    input  wire         start,
    output wire [N-1:0] out,
    output wire         nar,
    output wire         zero,
    output wire         done
);

    localparam integer FRACW  = N - es - 1;
    localparam integer SIGW   = N - es + 1;
    localparam integer WW     = 2 * SIGW;
    localparam integer SCALEW = Bs + es + 2;
    localparam integer LOWW   = WW - SIGW;

    localparam [SCALEW:0] WW_CONST = WW;
    localparam signed [SCALEW:0] FORMAT_MAX_SCALE =
        (rS * (1 << es)) - 1;
    localparam signed [SCALEW:0] FORMAT_MIN_SCALE =
        -(rS * (1 << es));

    localparam [SCALEW-1:0] SCALE_MAX_BITS =
        {1'b0, {(SCALEW-1){1'b1}}};
    localparam [SCALEW-1:0] SCALE_MIN_BITS =
        {1'b1, {(SCALEW-1){1'b0}}};

`ifndef SYNTHESIS
    initial begin
        if (N < 4)       $display("ERROR: N must be at least 4");
        if (es < 1)      $display("ERROR: es must be at least 1");
        if (rS < 2)      $display("ERROR: rS must be at least 2");
        if (rS > N-2)    $display("ERROR: rS must be <= N-2");
        if (N < rS+es+2) $display("ERROR: require N >= rS + es + 2");
        if ((1 << Bs) < N)
            $display("ERROR: Bs is too small for N");
    end
`endif

    assign done = start;

    wire sa = in_a[N-1];
    wire sb = in_b[N-1];

    wire body_nonzero_a = |in_a[N-2:0];
    wire body_nonzero_b = |in_b[N-2:0];

    wire nar_a  =  sa & ~body_nonzero_a;
    wire nar_b  =  sb & ~body_nonzero_b;
    wire zero_a = ~sa & ~body_nonzero_a;
    wire zero_b = ~sb & ~body_nonzero_b;

    wire any_nar = nar_a | nar_b;
    assign nar = any_nar;

    wire [N-1:0] mag_a = sa
        ? (~in_a + {{(N-1){1'b0}}, 1'b1})
        : in_a;

    wire [N-1:0] mag_b = sb
        ? (~in_b + {{(N-1){1'b0}}, 1'b1})
        : in_b;

    wire rc_a, rc_b;
    wire [Bs:0] regime_a, regime_b;
    wire [es-1:0] exp_a, exp_b;
    wire [FRACW-1:0] frac_a, frac_b;
    wire dec_zero_a, dec_zero_b;
    wire frac_zero_a, frac_zero_b;

    bposit_decode #(
        .N(N), .es(es), .rS(rS), .Bs(Bs)
    ) decode_a (
        .in(mag_a), .sign(), .is_zero(dec_zero_a), .is_nar(),
        .rc(rc_a), .regime(regime_a), .exp(exp_a),
        .frac(frac_a), .frac_zero(frac_zero_a)
    );

    bposit_decode #(
        .N(N), .es(es), .rS(rS), .Bs(Bs)
    ) decode_b (
        .in(mag_b), .sign(), .is_zero(dec_zero_b), .is_nar(),
        .rc(rc_b), .regime(regime_b), .exp(exp_b),
        .frac(frac_b), .frac_zero(frac_zero_b)
    );

    wire signed [Bs+1:0] regime_a_mag = {1'b0, regime_a};
    wire signed [Bs+1:0] regime_b_mag = {1'b0, regime_b};

    wire signed [Bs+1:0] regime_a_signed =
        rc_a ? regime_a_mag : -regime_a_mag;

    wire signed [Bs+1:0] regime_b_signed =
        rc_b ? regime_b_mag : -regime_b_mag;

    wire signed [SCALEW-1:0] scale_a = {regime_a_signed, exp_a};
    wire signed [SCALEW-1:0] scale_b = {regime_b_signed, exp_b};

    // Encoder-width-aligned significands.
    wire [SIGW-1:0] sig_a = {1'b1, frac_a, 1'b0};
    wire [SIGW-1:0] sig_b = {1'b1, frac_b, 1'b0};

    wire a_is_bigger =
          (scale_a > scale_b)
        | ((scale_a == scale_b) & (sig_a >= sig_b));

    wire result_sign_big   = a_is_bigger ? sa : sb;
    wire result_sign_small = a_is_bigger ? sb : sa;

    wire signed [SCALEW-1:0] scale_big =
        a_is_bigger ? scale_a : scale_b;

    wire signed [SCALEW-1:0] scale_small =
        a_is_bigger ? scale_b : scale_a;

    wire [SIGW-1:0] sig_big = a_is_bigger ? sig_a : sig_b;
    wire [SIGW-1:0] sig_small = a_is_bigger ? sig_b : sig_a;

    wire same_sign = (result_sign_big == result_sign_small);

    wire signed [SCALEW:0] scale_big_ext =
        {scale_big[SCALEW-1], scale_big};

    wire signed [SCALEW:0] scale_small_ext =
        {scale_small[SCALEW-1], scale_small};

    wire [SCALEW:0] scale_delta = scale_big_ext - scale_small_ext;
    wire shift_all = (scale_delta >= WW_CONST);
    wire [SCALEW:0] shift_amount = shift_all ? WW_CONST : scale_delta;

    wire [WW-1:0] sig_big_ext   = {sig_big,   {SIGW{1'b0}}};
    wire [WW-1:0] sig_small_ext = {sig_small, {SIGW{1'b0}}};

    wire [WW-1:0] shifted_small = shift_all
        ? {WW{1'b0}}
        : (sig_small_ext >> shift_amount);

    wire [WW-1:0] dropped_mask =
        (shift_amount == 0) ? {WW{1'b0}} :
        (shift_all)         ? {WW{1'b1}} :
                              ({WW{1'b1}} >> (WW - shift_amount));

    wire alignment_sticky = |(sig_small_ext & dropped_mask);

    wire [WW:0] add_full =
        {1'b0, sig_big_ext} + {1'b0, shifted_small};

    wire [WW-1:0] sub_raw = sig_big_ext - shifted_small;

    wire [WW-1:0] sub_corrected = alignment_sticky
        ? (sub_raw - {{(WW-1){1'b0}}, 1'b1})
        : sub_raw;

    wire add_carry = add_full[WW];

    wire [WW-1:0] add_normalized = add_carry
        ? add_full[WW:1]
        : add_full[WW-1:0];

    wire add_sticky = add_carry
        ? (add_full[0] | alignment_sticky)
        : alignment_sticky;

    wire signed [SCALEW:0] add_scale_ext =
        scale_big_ext + {{SCALEW{1'b0}}, add_carry};

    wire subtraction_is_zero =
        (sub_raw == {WW{1'b0}}) & ~alignment_sticky;

    reg [SCALEW:0] leading_zero_count;
    integer bit_index;

    always @(*) begin
        leading_zero_count = WW_CONST;
        for (bit_index = WW-1;
             bit_index >= 0;
             bit_index = bit_index - 1) begin
            if (sub_corrected[bit_index] &&
                (leading_zero_count == WW_CONST)) begin
                leading_zero_count = WW - 1 - bit_index;
            end
        end
    end

    wire [WW-1:0] sub_normalized = subtraction_is_zero
        ? {WW{1'b0}}
        : (sub_corrected << leading_zero_count);

    wire signed [SCALEW:0] leading_zero_signed = leading_zero_count;

    wire signed [SCALEW:0] sub_scale_ext =
        scale_big_ext - leading_zero_signed;

    wire [WW-1:0] result_sig_ext = same_sign
        ? add_normalized
        : sub_normalized;

    wire result_low_sticky = same_sign
        ? add_sticky
        : alignment_sticky;

    wire signed [SCALEW:0] result_scale_ext = same_sign
        ? add_scale_ext
        : sub_scale_ext;

    wire result_zero_math = ~same_sign & subtraction_is_zero;

    wire [SIGW-1:0] final_mantissa =
        result_sig_ext[WW-1 -: SIGW];

    wire final_sticky =
        result_low_sticky | (|result_sig_ext[LOWW-1:0]);

    wire signed [SCALEW:0] scale_max_ext =
        {1'b0, SCALE_MAX_BITS};

    wire signed [SCALEW:0] scale_min_ext =
        {1'b1, SCALE_MIN_BITS};

    wire signed [SCALEW-1:0] final_scale =
        (result_scale_ext > scale_max_ext) ? SCALE_MAX_BITS :
        (result_scale_ext < scale_min_ext) ? SCALE_MIN_BITS :
                                             result_scale_ext[SCALEW-1:0];

    wire [N-1:0] sum_encoded_raw;

    bposit_encode #(
        .N(N), .es(es), .rS(rS), .Bs(Bs)
    ) encode_sum (
        .sign(result_sign_big),
        .scale(final_scale),
        .mant_normalized(final_mantissa),
        .sticky_in(final_sticky),
        .special_nar(any_nar),
        .special_zero(result_zero_math),
        .out(sum_encoded_raw)
    );

    localparam [N-1:0] POS_MAX = {1'b0, {(N-1){1'b1}}};
    localparam [N-1:0] POS_MIN = {{(N-1){1'b0}}, 1'b1};

    wire [N-1:0] NEG_MAX =
        ~POS_MAX + {{(N-1){1'b0}}, 1'b1};

    wire [N-1:0] NEG_MIN =
        ~POS_MIN + {{(N-1){1'b0}}, 1'b1};

    wire [N-1:0] signed_max_result =
        result_sign_big ? NEG_MAX : POS_MAX;

    wire [N-1:0] signed_min_result =
        result_sign_big ? NEG_MIN : POS_MIN;

    wire result_overflow  = result_scale_ext > FORMAT_MAX_SCALE;
    wire result_underflow = result_scale_ext < FORMAT_MIN_SCALE;

    wire sum_collapsed_to_zero =
        (sum_encoded_raw == {N{1'b0}}) &
        ~result_zero_math &
        ~any_nar;

    wire [N-1:0] sum_encoded =
        result_zero_math       ? {N{1'b0}} :
        result_overflow        ? signed_max_result :
        result_underflow       ? signed_min_result :
        sum_collapsed_to_zero  ? signed_min_result :
                                 sum_encoded_raw;

    wire output_is_zero =
          (zero_a & zero_b)
        | (~zero_a & ~zero_b & result_zero_math);

    assign zero = output_is_zero & ~any_nar;

    assign out = any_nar ? {1'b1, {(N-1){1'b0}}} :
                 zero_a  ? in_b :
                 zero_b  ? in_a :
                           sum_encoded;

endmodule

// ============================================================================
// Sequential normal unfused EXACT bounded-posit MAC
//
// Stage 1: approximate multiplication -> product register
// Stage 2: bounded-posit addition      -> output register
//
// Latency: two rising clock edges from accepted input to valid output.
// Throughput: one new operation per clock after pipeline fill.
// Reset: active-high synchronous reset.
// ============================================================================
module bmac #(
    parameter integer N  = 12,
    parameter integer es = 2
)(
    input  wire         clk,
    input  wire         rst,
    input  wire         valid_in,

    input  wire [N-1:0] in1,
    input  wire [N-1:0] in2,
    input  wire [N-1:0] in3,

    output reg          valid_out,
    output reg  [N-1:0] out,
    output reg          nar,
    output reg          zero
);

    localparam integer rS = 4;
    localparam integer Bs =
        (N <= 8)  ? 3 :
        (N <= 16) ? 4 :
        (N <= 32) ? 5 :
        (N <= 64) ? 6 : 7;

    // Stage 1 combinational exact multiplier.
    wire [N-1:0] mult_result;

    posit_mult_exact #(
        .N(N),
        .es(es),
        .rS(rS),
        .Bs(Bs)
    ) stage1_multiplier (
        .in1(in1),
        .in2(in2),
        .start(valid_in),
        .out(mult_result),
        .nar(),
        .zero(),
        .done()
    );

    // Pipeline registers between multiplication and addition.
    reg [N-1:0] product_reg;
    reg [N-1:0] psum_reg;
    reg         valid_stage1;

    // Stage 2 combinational bounded-posit adder.
    wire [N-1:0] add_result;
    wire         add_nar;
    wire         add_zero;

    bposit_add #(
        .N(N),
        .es(es),
        .rS(rS),
        .Bs(Bs)
    ) stage2_adder (
        .in_a(product_reg),
        .in_b(psum_reg),
        .start(valid_stage1),
        .out(add_result),
        .nar(add_nar),
        .zero(add_zero),
        .done()
    );

    always @(posedge clk) begin
        if (rst) begin
            product_reg  <= {N{1'b0}};
            psum_reg     <= {N{1'b0}};
            valid_stage1 <= 1'b0;

            out       <= {N{1'b0}};
            nar       <= 1'b0;
            zero      <= 1'b1;
            valid_out <= 1'b0;
        end else begin
            // Stage 1: accept a new operation.
            valid_stage1 <= valid_in;
            if (valid_in) begin
                product_reg <= mult_result;
                psum_reg    <= in3;
            end

            // Stage 2: register the corresponding MAC result.
            valid_out <= valid_stage1;
            if (valid_stage1) begin
                out  <= add_result;
                nar  <= add_nar;
                zero <= add_zero;
            end
        end
    end

endmodule

module bsingalPE #(
    parameter integer N  = 16,
    parameter integer es = 2
)(
    input  wire         clk,
    input  wire         rst,

    input  wire         load_weight,
    input  wire [N-1:0] weight_in,

    input  wire         valid_in,
    input  wire [N-1:0] activation_in,
    input  wire [N-1:0] psum_in,

    output wire         valid_out,
    output wire [N-1:0] activation_out,
    output wire [N-1:0] psum_out,
    output wire         nar_out,
    output wire         zero_out
);

    reg [N-1:0] weight_reg;
    reg [N-1:0] activation_d1;
    reg [N-1:0] activation_d2;

    wire [N-1:0] mac_out;
    wire         mac_valid_out;
    wire         mac_nar;
    wire         mac_zero;

    bmac #(
        .N(N),
        .es(es)
    ) mac_core (
        .clk(clk),
        .rst(rst),
        .valid_in(valid_in),
        .in1(activation_in),
        .in2(weight_reg),
        .in3(psum_in),
        .valid_out(mac_valid_out),
        .out(mac_out),
        .nar(mac_nar),
        .zero(mac_zero)
    );

    assign valid_out      = mac_valid_out;
    assign activation_out = activation_d2;
    assign psum_out       = mac_out;
    assign nar_out        = mac_nar;
    assign zero_out       = mac_zero;

    always @(posedge clk) begin
        if (rst) begin
            weight_reg    <= {N{1'b0}};
            activation_d1 <= {N{1'b0}};
            activation_d2 <= {N{1'b0}};
        end else begin
            if (load_weight)
                weight_reg <= weight_in;

            activation_d1 <= activation_in;
            activation_d2 <= activation_d1;
        end
    end

`ifndef SYNTHESIS
    always @(posedge clk) begin
        if (!rst && load_weight && valid_in)
            $display("WARNING: load_weight and valid_in are both high; the current operation uses the previously stored weight.");
    end
`endif

endmodule
