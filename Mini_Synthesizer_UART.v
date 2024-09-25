`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: ECE 3300L
// Engineer: Mauricio Herrera
// 
// Create Date: 04/27/2024 06:27:54 PM
// Design Name: Lab 8
// Module Name: Mini_Synthesizer_UART
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module Mini_Synthesizer_UART(
   input clk,
   input rreset,
   output ampPWM,
   output ampSD,
   input UART_RX,
   input [15:0] SW,
   output [4:0] LED
    );
   
   wire pwm;
   wire	rx_dv;
   wire note_start;
   wire [7:0] sinpos1, envout;
   wire [7:0] rx_byte;
   wire [10:0] freq;

   assign ampPWM = (pwm) ? 1'bz : 1'b0;
   assign ampSD = 1'b1;

   button_pulse u1 (.clk(clk), .raw_button(rreset), .button_pulse(reset));
   uart_rx_vlog rx0 (clk, UART_RX, rx_dv, rx_byte);
   rx_parse_note rxp (clk, reset, rx_dv, rx_byte, freq, note_start);
   // this example sets the volume to drop from E0 to 00 in 200 cycles
   envel_gen egen (clk, reset, SW[15:12], SW[11:8], SW[7:4], SW[3:0], note_start, envout, LED);
   sine_gen sigen(clk, reset, (freq), envout, sinpos1);
   pwm_gen pwmg(clk, reset, sinpos1, pwm);
    
endmodule


module envel_gen(clk, reset, attack_time, decay_time, sustain_amp, release_time, start, envout, led);
   input clk, reset;
   input start;
   input [3:0] attack_time, decay_time, sustain_amp, release_time;
   output reg [7:0] envout;
   output reg [4:0] led;
   reg running;
   reg [7:0]  beginv, endv;
   reg [10:0] N;
   reg [10:0] Ncnt;
   reg [15:0] clkcnt;
   wire [31:0] m, mtop, mx, envout_calc;
 
   reg dstart, mstart1;
   wire dcomp, dovf, mcomp1, movf1;

// calculate the slope m = (endv - beginv) / N
   
   qadd #(20,32) subm ({1'b0, {3'b0,endv}, 20'b0}, {1'b1, {3'b0,beginv}, 20'b0}, mtop);
   qdiv #(20,32) mcalc (mtop, {1'b0, N, 20'b0}, dstart, clk, m, dcomp, dovf);

// calculate the value at cycle Ncnt: envout = m * Ncnt + beginv

   qmults #(20,32) mxcalc (m, {1'b0, Ncnt, 20'b0}, mstart1, clk, mx, mcomp1, movf1);
   qadd #(20,32) addb (mx, {1'b0, {3'b0,beginv}, 20'b0}, envout_calc);
   
   parameter IDLE = 0;
   parameter ATTACK = 1;
   parameter DECAY = 2;
   parameter SUSTAIN = 3;
   parameter RELEASE = 4;
   
   reg [2:0] state;
   reg [2:0] next_state;
   reg [7:0] next_endv, next_beginv;
   reg [10:0] next_N;
   reg next_running;
   
   always@(posedge clk or posedge reset)
     begin
        if (reset) begin
            state <= IDLE;
            envout <= 0;
            clkcnt <= 0;
            Ncnt <= 0;
            dstart <= 0;
            mstart1 <= 0;
        end else begin
            dstart <= 0;
            mstart1 <= 0;
            state <= next_state;
            endv <= next_endv;
            beginv <= next_beginv;
            N <= next_N;
            running <= next_running;
            if(start) begin
                Ncnt <= 0;
                clkcnt <= 0;
                //envout <= beginv;     
            end else if (running) begin
                clkcnt <= clkcnt + 1;
                if (clkcnt == 16'hFFFF) begin
                   if (Ncnt >= N) begin
                     // envout <= endv;
                      clkcnt <= 0;    
                      Ncnt <= 0;
                   end else begin
                      envout <= envout_calc[27:20];
                      Ncnt <= Ncnt + 1;
                   end
                end
                if (clkcnt == 16'h8100) begin
                   mstart1 <= 1;  // recalculte output
                end
                if (Ncnt == 0 && clkcnt == 16'h0001) begin 
                   dstart <= 1;   // recalculate m once per envelope start
                end
            end 
	     end
     end
     
     always @(*) begin
     next_state = state;
     next_endv = endv;
     next_beginv = beginv;
     next_N = N;
     next_running = running;
     if (reset)
        next_state = IDLE;
     else begin
     case(state)
        IDLE: begin
            led = 5'b00000;
            next_running = 0;
            next_beginv = 0;
            next_endv = 0;
            next_N = 0;
            if (start) begin
                next_running = 1;
                next_state = ATTACK;
            end
        end
        ATTACK: begin
        if(start)
        next_state = ATTACK;
        else begin
            led = 5'b00001;
            next_N = {attack_time, 7'b0000000};
            next_beginv = 8'h00;
            next_endv = 8'hE0;
                if (clkcnt == 16'hFFFF && Ncnt >= next_N)begin
                    next_state = DECAY;  
                    end
        end
        end
        DECAY: begin
        if(start)
        next_state = ATTACK;
        else begin
            led = 5'b00010;
            next_N = {decay_time, 7'b0000000};
            next_beginv = 8'hE0;
            next_endv = {sustain_amp, 4'b0000};
                if (clkcnt == 16'hFFFF && Ncnt >= next_N)begin
                    next_state = SUSTAIN;
                end
        end
        end
        SUSTAIN: begin
        if(start) 
            next_state = ATTACK;
        else begin
            led = 5'b00100;
            next_N = {11'b10111110101};
            next_beginv = {sustain_amp, 4'b0000};
            next_endv = {sustain_amp, 4'b0000};
                if (clkcnt == 16'hFFFF && Ncnt >= next_N)begin
                    next_state = RELEASE;
                end
        end
        end
        RELEASE: begin
        if(start)
            next_state = ATTACK;
        else begin
            led = 5'b01000;
            next_N = {release_time, 7'b0000000};
            next_beginv = {sustain_amp, 4'b0000};
            next_endv = 8'h00;
                if (clkcnt == 16'hFFFF && Ncnt >= next_N)begin
                     next_running = 0;
                     next_state = IDLE; 
                end
        end
        end
        default: begin
        end
     endcase
     
     end
     end
endmodule


module pwm_gen(clk, reset, inp, pwm); 
   input clk, reset;
   input [7:0] inp;
   output pwm;

   reg [7:0] pwmcnt;
   wire	       pwm;

   assign pwm = (pwmcnt < inp) ? 1 : 0;

  always@(posedge clk or posedge reset)
     begin
         if (reset) begin
	        pwmcnt <= 0;
         end else begin
	        pwmcnt <= pwmcnt + 1;
	     end
     end

endmodule

module rx_parse_note(
   input clk,
   input reset, 
   input rx_dv,
   input [7:0] rx_byte, 
   output reg [10:0] freq,
   output reg note_start
   );

    always @(posedge clk or posedge reset)
    begin
        if (reset) begin
           freq <= 888;
           note_start <= 0;
        end else begin
           if (rx_dv) begin
             note_start <= 1;
             case (rx_byte)
               "a": freq <= 1493; // 261.6Hz (C4)
               "w": freq <= 1409; // 277.2Hz (C#/Db)
               "s": freq <= 1330; // 293.7Hz (D)
               "e": freq <= 1256; // 311.1Hz (D#/Eb)
               "d": freq <= 1185; // 329.6Hz (E)
               "f": freq <= 1119; // 349.2Hz (F)
               "t": freq <= 1056; // 370Hz   (F#/Gb)
               "g": freq <= 996;  // 392Hz   (G)
               "y": freq <= 941;  // 415.3Hz (G#/Ab)
               "h": freq <= 888;  // 440Hz   (A)
               "u": freq <= 838;  // 466.2Hz (A#/Bb)
               "j": freq <= 791;  // 493.9Hz (B)
               "k": freq <= 746;  // 523.2Hz (C5)
               "o": freq <= 704;  // 554.4Hz (C#/Db5)
               "l": freq <= 665;  // 587.3Hz (D5)
               "p": freq <= 628;  // 622.3Hz (D#/Eb5)
               ";": freq <= 593;  // 659.3Hz (E5)
               "'": freq <= 559;  // 698.5Hz (F5)
               "[": freq <= 528;  // 739.9Hz (F#/Gb5)
               "]": freq <= 498;  // 783.9Hz (G5)
             endcase 
           end else begin
             note_start <= 0;
           end
        end
    end

endmodule

module uart_rx_vlog 
  (
   input        i_Clock,           // master clock
   input        i_Rx_Serial,       // UART serial RX wire
   output       o_Rx_DV,           // data valid - goes high for 1 clk cycle when data byte has been received
   output [7:0] o_Rx_Byte          // byte received
   );

  localparam CLKS_PER_BIT = 14'b10100010110000;   
  localparam s_IDLE         = 3'b000;
  localparam s_RX_START_BIT = 3'b001;
  localparam s_RX_DATA_BITS = 3'b010;
  localparam s_RX_STOP_BIT  = 3'b011;
  localparam s_CLEANUP      = 3'b100;
   
  reg           r_Rx_Data_R = 1'b1;
  reg           r_Rx_Data   = 1'b1;
   
  reg [13:0]     r_Clock_Count = 0;
  reg [2:0]     r_Bit_Index   = 0; //8 bits total
  reg [7:0]     r_Rx_Byte     = 0;
  reg           r_Rx_DV       = 0;
  reg [2:0]     r_SM_Main     = 0;
   
  // Purpose: Double-register the incoming data.
  // This allows it to be used in the UART RX Clock Domain.
  // (It removes problems caused by metastability)
  always @(posedge i_Clock)
    begin
      r_Rx_Data_R <= i_Rx_Serial;
      r_Rx_Data   <= r_Rx_Data_R;
    end
   
   
  // Purpose: Control RX state machine
  always @(posedge i_Clock)
    begin
       
      case (r_SM_Main)
        s_IDLE :
          begin
            r_Rx_DV       <= 1'b0;
            r_Clock_Count <= 0;
            r_Bit_Index   <= 0;
             
            if (r_Rx_Data == 1'b0)          // Start bit detected
              r_SM_Main <= s_RX_START_BIT;
            else
              r_SM_Main <= s_IDLE;
          end
         
        // Check middle of start bit to make sure it's still low
        s_RX_START_BIT :
          begin
            if (r_Clock_Count == (CLKS_PER_BIT-1)/2)
              begin
                if (r_Rx_Data == 1'b0)
                  begin
                    r_Clock_Count <= 0;  // reset counter, found the middle
                    r_SM_Main     <= s_RX_DATA_BITS;
                  end
                else
                  r_SM_Main <= s_IDLE;
              end
            else
              begin
                r_Clock_Count <= r_Clock_Count + 1;
                r_SM_Main     <= s_RX_START_BIT;
              end
          end // case: s_RX_START_BIT
         
         
        // Wait CLKS_PER_BIT-1 clock cycles to sample serial data
        s_RX_DATA_BITS :
          begin
            if (r_Clock_Count < CLKS_PER_BIT-1)
              begin
                r_Clock_Count <= r_Clock_Count + 1;
                r_SM_Main     <= s_RX_DATA_BITS;
              end
            else
              begin
                r_Clock_Count          <= 0;
                r_Rx_Byte[r_Bit_Index] <= r_Rx_Data;
                 
                // Check if we have received all bits
                if (r_Bit_Index < 7)
                  begin
                    r_Bit_Index <= r_Bit_Index + 1;
                    r_SM_Main   <= s_RX_DATA_BITS;
                  end
                else
                  begin
                    r_Bit_Index <= 0;
                    r_SM_Main   <= s_RX_STOP_BIT;
                  end
              end
          end // case: s_RX_DATA_BITS
     
     
        // Receive Stop bit.  Stop bit = 1
        s_RX_STOP_BIT :
          begin
            // Wait CLKS_PER_BIT-1 clock cycles for Stop bit to finish
            if (r_Clock_Count < CLKS_PER_BIT-1)
              begin
                r_Clock_Count <= r_Clock_Count + 1;
                r_SM_Main     <= s_RX_STOP_BIT;
              end
            else
              begin
                r_Rx_DV       <= 1'b1;
                r_Clock_Count <= 0;
                r_SM_Main     <= s_CLEANUP;
              end
          end // case: s_RX_STOP_BIT
     
         
        // Stay here 1 clock
        s_CLEANUP :
          begin
            r_SM_Main <= s_IDLE;
            r_Rx_DV   <= 1'b0;
          end
         
         
        default :
          r_SM_Main <= s_IDLE;
         
      endcase
    end   
   
  assign o_Rx_DV   = r_Rx_DV;
  assign o_Rx_Byte = r_Rx_Byte;
   
endmodule // uart_rx

module button_pulse(
		    input clk,
		    input raw_button,
		    output button_pulse
		    );
    
    localparam N = 3;
    
    reg [N - 1:0] Q_reg;
    
    always @(posedge clk)
    begin
        Q_reg <= {Q_reg[N - 2:0], raw_button};
    end
    
    assign button_pulse = (&Q_reg[N - 2:0]) & ~Q_reg[N-1];
endmodule

module qadd #(
	//Parameterized values
	parameter Q = 15,
	parameter N = 32
	)
	(
    input [N-1:0] a,
    input [N-1:0] b,
    output [N-1:0] c
    );

reg [N-1:0] res;

assign c = res;

always @(a,b) begin
	// both negative or both positive
	if(a[N-1] == b[N-1]) begin						//	Since they have the same sign, absolute magnitude increases
		res[N-2:0] = a[N-2:0] + b[N-2:0];		//		So we just add the two numbers
		res[N-1] = a[N-1];							//		and set the sign appropriately...  Doesn't matter which one we use, 
															//		they both have the same sign
															//	Do the sign last, on the off-chance there was an overflow...  
		end												//		Not doing any error checking on this...
	//	one of them is negative...
	else if(a[N-1] == 0 && b[N-1] == 1) begin		//	subtract a-b
		if( a[N-2:0] > b[N-2:0] ) begin					//	if a is greater than b,
			res[N-2:0] = a[N-2:0] - b[N-2:0];			//		then just subtract b from a
			res[N-1] = 0;										//		and manually set the sign to positive
			end
		else begin												//	if a is less than b,
			res[N-2:0] = b[N-2:0] - a[N-2:0];			//		we'll actually subtract a from b to avoid a 2's complement answer
			if (res[N-2:0] == 0)
				res[N-1] = 0;										//		I don't like negative zero....
			else
				res[N-1] = 1;										//		and manually set the sign to negative
			end
		end
	else begin												//	subtract b-a (a negative, b positive)
		if( a[N-2:0] > b[N-2:0] ) begin					//	if a is greater than b,
			res[N-2:0] = a[N-2:0] - b[N-2:0];			//		we'll actually subtract b from a to avoid a 2's complement answer
			if (res[N-2:0] == 0)
				res[N-1] = 0;										//		I don't like negative zero....
			else
				res[N-1] = 1;										//		and manually set the sign to negative
			end
		else begin												//	if a is less than b,
			res[N-2:0] = b[N-2:0] - a[N-2:0];			//		then just subtract a from b
			res[N-1] = 0;										//		and manually set the sign to positive
			end
		end
	end
endmodule

module qdiv #(
	//Parameterized values
	parameter Q = 15,
	parameter N = 32
	)
	(
	input 	[N-1:0] i_dividend,
	input 	[N-1:0] i_divisor,
	input 	i_start,
	input 	i_clk,
	output 	[N-1:0] o_quotient_out,
	output 	o_complete,
	output	o_overflow
	);
 
	reg [2*N+Q-3:0]	reg_working_quotient;	//	Our working copy of the quotient
	reg [N-1:0] 		reg_quotient;				//	Final quotient
	reg [N-2+Q:0] 		reg_working_dividend;	//	Working copy of the dividend
	reg [2*N+Q-3:0]	reg_working_divisor;		// Working copy of the divisor
 
	reg [N-1:0] 			reg_count; 		//	This is obviously a lot bigger than it needs to be, as we only need 
													//		count to N-1+Q but, computing that number of bits requires a 
													//		logarithm (base 2), and I don't know how to do that in a 
													//		way that will work for everyone
										 
	reg					reg_done;			//	Computation completed flag
	reg					reg_sign;			//	The quotient's sign bit
	reg					reg_overflow;		//	Overflow flag
 
	initial reg_done = 1'b1;				//	Initial state is to not be doing anything
	initial reg_overflow = 1'b0;			//		And there should be no woverflow present
	initial reg_sign = 1'b0;				//		And the sign should be positive

	initial reg_working_quotient = 0;	
	initial reg_quotient = 0;				
	initial reg_working_dividend = 0;	
	initial reg_working_divisor = 0;		
 	initial reg_count = 0; 		

 
	assign o_quotient_out[N-2:0] = reg_quotient[N-2:0];	//	The division results
	assign o_quotient_out[N-1] = reg_sign;						//	The sign of the quotient
	assign o_complete = reg_done;
	assign o_overflow = reg_overflow;
 
	always @( posedge i_clk ) begin
		if( reg_done && i_start ) begin										//	This is our startup condition
			//  Need to check for a divide by zero right here, I think....
			reg_done <= 1'b0;												//	We're not done			
			reg_count <= N+Q-1;											//	Set the count
			reg_working_quotient <= 0;									//	Clear out the quotient register
			reg_working_dividend <= 0;									//	Clear out the dividend register 
			reg_working_divisor <= 0;									//	Clear out the divisor register 
			reg_overflow <= 1'b0;										//	Clear the overflow register

			reg_working_dividend[N+Q-2:Q] <= i_dividend[N-2:0];				//	Left-align the dividend in its working register
			reg_working_divisor[2*N+Q-3:N+Q-1] <= i_divisor[N-2:0];		//	Left-align the divisor into its working register

			reg_sign <= i_dividend[N-1] ^ i_divisor[N-1];		//	Set the sign bit
			end 
		else if(!reg_done) begin
			reg_working_divisor <= reg_working_divisor >> 1;	//	Right shift the divisor (that is, divide it by two - aka reduce the divisor)
			reg_count <= reg_count - 1;								//	Decrement the count

			//	If the dividend is greater than the divisor
			if(reg_working_dividend >= reg_working_divisor) begin
				reg_working_quotient[reg_count] <= 1'b1;										//	Set the quotient bit
				reg_working_dividend <= reg_working_dividend - reg_working_divisor;	//		and subtract the divisor from the dividend
				end
 
			//stop condition
			if(reg_count == 0) begin
				reg_done <= 1'b1;										//	If we're done, it's time to tell the calling process
				reg_quotient <= reg_working_quotient;			//	Move in our working copy to the outside world
				if (reg_working_quotient[2*N+Q-3:N]>0)
					reg_overflow <= 1'b1;
					end
			else
				reg_count <= reg_count - 1;	
			end
		end
endmodule

module qmult #(
	//Parameterized values
	parameter Q = 15,
	parameter N = 32
	)
	(
	 input			[N-1:0]	i_multiplicand,
	 input			[N-1:0]	i_multiplier,
	 output			[N-1:0]	o_result,
	 output	reg				ovr
	 );
	 
	 //	The underlying assumption, here, is that both fixed-point values are of the same length (N,Q)
	 //		Because of this, the results will be of length N+N = 2N bits....
	 //		This also simplifies the hand-back of results, as the binimal point 
	 //		will always be in the same location...
	
	reg [2*N-1:0]	r_result;		//	Multiplication by 2 values of N bits requires a 
											//		register that is N+N = 2N deep...
	reg [N-1:0]		r_RetVal;
	
//--------------------------------------------------------------------------------
	assign o_result = r_RetVal;	//	Only handing back the same number of bits as we received...
											//		with fixed point in same location...
	
//---------------------------------------------------------------------------------
	always @(i_multiplicand, i_multiplier)	begin						//	Do the multiply any time the inputs change
		r_result <= i_multiplicand[N-2:0] * i_multiplier[N-2:0];	//	Removing the sign bits from the multiply - that 
																					//		would introduce *big* errors	
		ovr <= 1'b0;															//	reset overflow flag to zero
		end
	
		//	This always block will throw a warning, as it uses a & b, but only acts on changes in result...
	always @(r_result) begin													//	Any time the result changes, we need to recompute the sign bit,
		r_RetVal[N-1] <= i_multiplicand[N-1] ^ i_multiplier[N-1];	//		which is the XOR of the input sign bits...  (you do the truth table...)
		r_RetVal[N-2:0] <= r_result[N-2+Q:Q];								//	And we also need to push the proper N bits of result up to 
																						//		the calling entity...
		if (r_result[2*N-2:N-1+Q] > 0)										// And finally, we need to check for an overflow
			ovr <= 1'b1;
		end

endmodule

module qmults#(
	//Parameterized values
	parameter Q = 15,
	parameter N = 32
	)
	(
	input 	[N-1:0]  i_multiplicand,
	input 	[N-1:0]	i_multiplier,
	input 	i_start,
	input 	i_clk,
	output 	[N-1:0] o_result_out,
	output 	o_complete,
	output	o_overflow
	);

	reg [2*N-2:0]	reg_working_result;		//	a place to accumulate our result
	reg [2*N-2:0]	reg_multiplier_temp;		//	a working copy of the multiplier
	reg [N-1:0]		reg_multiplicand_temp;	//	a working copy of the umultiplicand
	
	reg [N-1:0] 			reg_count; 		//	This is obviously a lot bigger than it needs to be, as we only need 
												//		count to N, but computing that number of bits requires a 
												//		logarithm (base 2), and I don't know how to do that in a 
												//		way that will work for every possibility
										 
	reg					reg_done;		//	Computation completed flag
	reg					reg_sign;		//	The result's sign bit
	reg					reg_overflow;	//	Overflow flag
 
	initial reg_done = 1'b1;			//	Initial state is to not be doing anything
	initial reg_overflow = 1'b0;		//		And there should be no woverflow present
	initial reg_sign = 1'b0;			//		And the sign should be positive
	
	assign o_result_out[N-2:0] = reg_working_result[N-2+Q:Q];	//	The multiplication results
	assign o_result_out[N-1] = reg_sign;								//	The sign of the result
	assign o_complete = reg_done;											//	"Done" flag
	assign o_overflow = reg_overflow;									//	Overflow flag
	
	always @( posedge i_clk ) begin
		if( reg_done && i_start ) begin										//	This is our startup condition
			reg_done <= 1'b0;														//	We're not done			
			reg_count <= 0;														//	Reset the count
			reg_working_result <= 0;											//	Clear out the result register
			reg_multiplier_temp <= 0;											//	Clear out the multiplier register 
			reg_multiplicand_temp <= 0;										//	Clear out the multiplicand register 
			reg_overflow <= 1'b0;												//	Clear the overflow register

			reg_multiplicand_temp <= i_multiplicand[N-2:0];				//	Load the multiplicand in its working register and lose the sign bit
			reg_multiplier_temp <= i_multiplier[N-2:0];					//	Load the multiplier into its working register and lose the sign bit

			reg_sign <= i_multiplicand[N-1] ^ i_multiplier[N-1];		//	Set the sign bit
			end 

		else if (!reg_done) begin
			if (reg_multiplicand_temp[reg_count] == 1'b1)								//	if the appropriate multiplicand bit is 1
				reg_working_result <= reg_working_result + reg_multiplier_temp;	//		then add the temp multiplier
	
			reg_multiplier_temp <= reg_multiplier_temp << 1;						//	Do a left-shift on the multiplier
			reg_count <= reg_count + 1;													//	Increment the count

			//stop condition
			if(reg_count == N) begin
				reg_done <= 1'b1;										//	If we're done, it's time to tell the calling process
				if (reg_working_result[2*N-2:N-1+Q] > 0)			// Check for an overflow
					reg_overflow <= 1'b1;
//			else
//				reg_count <= reg_count + 1;													//	Increment the count
				end
			end
		end
endmodule

module sine_gen(clk, reset, N, ampl, sinpos);
   input clk, reset;
   input [10:0] N;  
   input [7:0] ampl;
   output [7:0] sinpos;

   wire [31:0] sin_out, cos_out;
   reg [31:0] sin_r, cos_r;
   reg [10:0] Ncnt;
   reg [7:0] pwmcnt;
   wire [31:0] a, delsin, delcos, sinfpos;
   wire [7:0] sinpos;

   reg	       dstart, mstart1, mstart2;
   wire	       dcomp, dovf, mcomp1, movf1, mcomp2, movf2;
   wire [6:0]  ampl_div_2;

   assign ampl_div_2 = ampl >> 1;  // the sine amplitude is 1/2 the peak-to-peak value

// a = (2 * Pi) / N
//   2 * Pi shifted left by 20 bits = 6588397 = 0x006_487ED in 20,32 format
//   N in 20,32 format: {1'b0, N, 20'b0}   

   qdiv #(20,32) acalc (32'h006_487ED, {1'b0, N, 20'b0}, dstart, clk, a, dcomp, dovf);

//   sin_out = sin_r + a*(cos_r);
//   cos_out = cos_r - a*(sin_r);
   
// its a bit wasteful to use two multipliers, we could multiplex into one,
// but we have plenty of FPGA resouces on this chip...   

   qmults #(20,32) dscalc (a, cos_r, mstart1, clk, delsin, mcomp1, movf1);

   qmults #(20,32) dccalc (a, sin_r, mstart2, clk, delcos, mcomp2, movf2);

   qadd #(20,32) adds (sin_r, delsin, sin_out);

   qadd #(20,32) subc (cos_r, {~delcos[31],delcos[30:0]}, cos_out);

// sin function is positive and negative
// add an offset so the sinpos is always positive
   
   qadd #(20,32) addos (sin_r, {1'b0, {4'b0,ampl_div_2} + 11'h1, 20'b0}, sinfpos);
   
   assign sinpos = sinfpos[31] ? 0 : sinfpos[27:20];  // output the integer part, clip to zero if negative
   
   always@(posedge clk or posedge reset)
     begin
         if (reset) begin
            sin_r <= 32'h000_00000;
            cos_r <= {1'b0, {4'b0,ampl_div_2}, 20'b0};
	        Ncnt <= 0;
	        pwmcnt <= 0;
	        dstart <= 0;
	        mstart1 <= 0;
	        mstart2 <= 0;
         end else begin
	        dstart <= 0;
	        mstart1 <= 0;
	        mstart2 <= 0;
	        pwmcnt <= pwmcnt + 1;

	      if (pwmcnt == 8'hFF) begin
	       if (Ncnt >= N-1) begin     
		      sin_r <= 32'h000_00000; // reset sin & cos for the first cycle,
		      cos_r <= {1'b0, {4'b0,ampl_div_2}, 20'b0}; //   so we don't accumulate errors
		      Ncnt <= 0;
	       end else begin
		      sin_r <= sin_out;
		      cos_r <= cos_out;
		      Ncnt <= Ncnt + 1;
	       end
	      end

	    if (pwmcnt == 8'h81) begin
	       mstart1 <= 1;  // recalculte cos & sin
	       mstart2 <= 1;
	    end

 	    if (Ncnt == 0 && pwmcnt == 8'h01) begin 
	       dstart <= 1;   // recalculate a once per audio cycle
	    end
	 end 
     end
endmodule