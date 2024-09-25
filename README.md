# FPGA-Mini-Synthesizer
For this project, I implemented a digital synthesizer on an FPGA using Verilog. The synthesizer plays a complete two-octave chromatic scale based on ASCII character input via UART and features amplitude modulation using an ADSR (Attack, Decay, Sustain, Release) envelope. This project features Sine Wave generation, Pulse Width Modulation, ADSR Envelope Modulation, Real-Time Control, and UART input for note control. LEDs are used for debugging and indicate the current state of the ADSR envelope when a new note is played on the keyboard.

The ADSR settings can be adjusted using the FPGA switches (SW 0-15).
SW[15:12]: Controls the attack time.
SW[11:8]: Sets the decay time.
SW[7:4]: Adjusts the sustain amplitude.
SW[3:0]: Sets the release time.
