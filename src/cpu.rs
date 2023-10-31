struct CpuFlags {
    pub carry: bool,
    pub zero: bool,
    pub int_dis: bool,
    pub dec_mode: bool,
    pub brk: bool,
    pub overflow: bool,
    pub negative: bool,
}

pub struct Cpu6502 {
    pc: u16,
    sp: u8,
    acc: u8,
    x: u8,
    y: u8,
    flags: CpuFlags,

    pub mem: [u8; 4096],
}

#[derive(Clone)]
enum AddressMode {
    Implicit,
    Accumulator,
    Immediate,
    ZeroPage,
    ZeroPageX,
    ZeroPageY,
    Relative,
    Absolute,
    AbsoluteX,
    AbsoluteY,
    Indirect,
    IndexedIndirect,
    IndirectIndexed,

    NA,
}

impl CpuFlags {
    fn new() -> Self {
        CpuFlags {
            carry: false,
            zero: false,
            int_dis: false,
            dec_mode: false,
            brk: false,
            overflow: false,
            negative: false,
        }
    }
}

impl Cpu6502 {
    pub fn new() -> Self {
        Cpu6502 {
            pc: 0,
            sp: 0,
            acc: 0,
            x: 0,
            y: 0,
            flags: CpuFlags::new(),
            mem: [0; 4096],
        }
    }

    pub fn execute(&mut self) -> Result<(), String> {
        let inst = self.read_mem_u8(self.pc);
        self.pc += 1;

        let a = ((inst >> 5) & 0x7) as usize;
        let b = ((inst >> 2) & 0x7) as usize;
        let c = (inst & 0x3) as usize;

        let (func, mode) = &DISPATCH_TABLE[c][a][b];
        func(self, mode.clone())?;

        Ok(())
    }

    // Add with Carry
    fn ex_adc(&mut self, mode: AddressMode) -> Result<(), String> {
        let mut result: u16 = (self.acc as u16) + self.get_operand(mode)?;
        if self.flags.carry {
            result += 1;
        }

        if self.acc == 0 {
            self.flags.zero = true;
        } else if (self.acc & 0x80) != 0 {
            self.flags.negative = true;
        }

        //TODO: don't know if this is right and there's probably a better way to do it
        if (result & 0xF0) != ((self.acc as u16) & 0xF0) {
            self.flags.overflow = true;
        }
        if (result & 0xFF00) != 0 {
            self.flags.carry = true;
        }

        self.acc = (result & 0xFF) as u8;
        Ok(())
    }

    fn ex_and(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_asl(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_bcc(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_bcs(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_beq(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_bit(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_bmi(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_bne(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_bpl(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_brk(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_bvc(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_bvs(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_clc(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_cld(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_cli(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_clv(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_cmp(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_cpx(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_cpy(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_dec(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_dex(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_dey(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_eor(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_inc(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_inx(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_iny(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_jmp(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_jsr(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_lda(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_ldx(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_ldy(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_lsr(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_nop(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_ora(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_pha(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_php(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_pla(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_plp(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_rol(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_ror(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_rti(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_rts(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_sbc(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_sec(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_sed(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_sei(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_sta(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_stx(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_sty(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_tax(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_tay(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_tsx(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_txa(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_txs(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    fn ex_tya(&mut self, _mode: AddressMode) -> Result<(), String> {
        Ok(())
    }

    // Unrecognized Opcode
    fn op_err(&mut self, _mode: AddressMode) -> Result<(), String> {
        Err(format!("Invalid Instruction"))
    }

    fn read_mem_u8(&self, addr: u16) -> u8 {
        self.mem[addr as usize]
    }

    fn read_mem_u16(&self, addr: u16) -> u16 {
        //little endian
        let mut result: u16 = (self.mem[(addr as usize) + 1] as u16) << 8;
        result |= self.mem[addr as usize] as u16;
        result
    }

    fn get_operand(&mut self, mode: AddressMode) -> Result<u16, String> {
        let m: u16;
        match mode {
            AddressMode::Immediate => {
                m = self.read_mem_u8(self.pc) as u16;
                self.pc += 1;
            }
            AddressMode::ZeroPage => {
                let temp = self.read_mem_u8(self.pc) as u16;
                self.pc += 1;
                m = self.read_mem_u8(temp) as u16;
            }
            AddressMode::ZeroPageX => {
                let temp = ((self.read_mem_u8(self.pc) as u16) + (self.x as u16)) & 0xFF;
                self.pc += 1;
                m = self.read_mem_u8(temp) as u16;
            }
            AddressMode::Absolute => {
                let temp = self.read_mem_u16(self.pc);
                self.pc += 2;
                m = self.read_mem_u8(temp) as u16;
            }
            AddressMode::AbsoluteX => {
                let temp = self.read_mem_u16(self.pc) + (self.x as u16);
                self.pc += 2;
                m = self.read_mem_u8(temp) as u16;
            }
            AddressMode::AbsoluteY => {
                let temp = self.read_mem_u16(self.pc) + (self.y as u16);
                self.pc += 2;
                m = self.read_mem_u8(temp) as u16;
            }
            AddressMode::IndexedIndirect => {
                let mut temp = (self.read_mem_u8(self.pc) as u16) + (self.x as u16);
                self.pc += 1;
                temp = self.read_mem_u16(temp & 0xFF);
                m = self.read_mem_u8(temp) as u16;
            }
            AddressMode::IndirectIndexed => {
                let mut temp = self.read_mem_u8(self.pc) as u16;
                self.pc += 1;
                temp = self.read_mem_u16(temp);
                temp += self.y as u16;
                m = self.read_mem_u8(temp) as u16;
            }
            _ => { return Err("Invalid address mode".to_string()); }
        }
        Ok(m)
    }
}

// Instruction Layout Reference: https://www.masswerk.at/6502/6502_instruction_set.html#layout
// Indexed by DISPATCH_TABLE[c][a][b]
static DISPATCH_TABLE: [[[(fn(&mut Cpu6502, AddressMode) -> Result<(), String>, AddressMode); 8]; 8]; 3] = [
    // c = 0
    [          // b = 0                                        = 1                                       = 2                                       = 3                                       = 4                                       = 5                                        = 6                                       = 7
    /* a = 0*/ [ (Cpu6502::ex_brk, AddressMode::Implicit),  (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_php, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_bpl, AddressMode::Relative), (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_clc, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::NA) ],
    /*     1*/ [ (Cpu6502::ex_jsr, AddressMode::Absolute),  (Cpu6502::ex_bit, AddressMode::ZeroPage), (Cpu6502::ex_plp, AddressMode::Implicit), (Cpu6502::ex_bit, AddressMode::Absolute), (Cpu6502::ex_bmi, AddressMode::Relative), (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_sec, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::NA) ],
    /*     2*/ [ (Cpu6502::ex_rti, AddressMode::Implicit),  (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_pha, AddressMode::Implicit), (Cpu6502::ex_jmp, AddressMode::Absolute), (Cpu6502::ex_bvc, AddressMode::Relative), (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_cli, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::NA) ],
    /*     3*/ [ (Cpu6502::ex_rts, AddressMode::Implicit),  (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_pla, AddressMode::Implicit), (Cpu6502::ex_jmp, AddressMode::Indirect), (Cpu6502::ex_bvs, AddressMode::Relative), (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_sei, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::NA) ],
    /*     4*/ [ (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_sty, AddressMode::ZeroPage), (Cpu6502::ex_dey, AddressMode::Implicit), (Cpu6502::ex_sty, AddressMode::Absolute), (Cpu6502::ex_bcc, AddressMode::Relative), (Cpu6502::ex_sty, AddressMode::ZeroPageX), (Cpu6502::ex_tya, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::NA) ],
    /*     5*/ [ (Cpu6502::ex_ldy, AddressMode::Immediate), (Cpu6502::ex_ldy, AddressMode::ZeroPage), (Cpu6502::ex_tay, AddressMode::Implicit), (Cpu6502::ex_ldy, AddressMode::Absolute), (Cpu6502::ex_bcs, AddressMode::Relative), (Cpu6502::ex_ldy, AddressMode::ZeroPageX), (Cpu6502::ex_clv, AddressMode::Implicit), (Cpu6502::ex_ldy, AddressMode::AbsoluteX) ],
    /*     6*/ [ (Cpu6502::ex_cpy, AddressMode::Immediate), (Cpu6502::ex_cpy, AddressMode::ZeroPage), (Cpu6502::ex_iny, AddressMode::Implicit), (Cpu6502::ex_cpy, AddressMode::Absolute), (Cpu6502::ex_bne, AddressMode::Relative), (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_cld, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::NA) ],
    /*     7*/ [ (Cpu6502::ex_cpx, AddressMode::Immediate), (Cpu6502::ex_cpx, AddressMode::ZeroPage), (Cpu6502::ex_inx, AddressMode::Implicit), (Cpu6502::ex_cpx, AddressMode::Absolute), (Cpu6502::ex_beq, AddressMode::Relative), (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_sed, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::NA) ],
    ],

    // c = 1
    [          // b = 0                                              = 1                                       = 2                                        = 3                                       = 4                                              = 5                                        = 6                                        = 7
    /* a = 0*/ [ (Cpu6502::ex_ora, AddressMode::IndexedIndirect), (Cpu6502::ex_ora, AddressMode::ZeroPage), (Cpu6502::ex_ora, AddressMode::Immediate), (Cpu6502::ex_ora, AddressMode::Absolute), (Cpu6502::ex_ora, AddressMode::IndirectIndexed), (Cpu6502::ex_ora, AddressMode::ZeroPageX), (Cpu6502::ex_ora, AddressMode::AbsoluteY), (Cpu6502::ex_ora, AddressMode::AbsoluteX) ],
    /*     1*/ [ (Cpu6502::ex_and, AddressMode::IndexedIndirect), (Cpu6502::ex_and, AddressMode::ZeroPage), (Cpu6502::ex_and, AddressMode::Immediate), (Cpu6502::ex_and, AddressMode::Absolute), (Cpu6502::ex_and, AddressMode::IndirectIndexed), (Cpu6502::ex_and, AddressMode::ZeroPageX), (Cpu6502::ex_and, AddressMode::AbsoluteY), (Cpu6502::ex_and, AddressMode::AbsoluteX) ],
    /*     2*/ [ (Cpu6502::ex_eor, AddressMode::IndexedIndirect), (Cpu6502::ex_eor, AddressMode::ZeroPage), (Cpu6502::ex_eor, AddressMode::Immediate), (Cpu6502::ex_eor, AddressMode::Absolute), (Cpu6502::ex_eor, AddressMode::IndirectIndexed), (Cpu6502::ex_eor, AddressMode::ZeroPageX), (Cpu6502::ex_eor, AddressMode::AbsoluteY), (Cpu6502::ex_eor, AddressMode::AbsoluteX) ],
    /*     3*/ [ (Cpu6502::ex_adc, AddressMode::IndexedIndirect), (Cpu6502::ex_adc, AddressMode::ZeroPage), (Cpu6502::ex_adc, AddressMode::Immediate), (Cpu6502::ex_adc, AddressMode::Absolute), (Cpu6502::ex_adc, AddressMode::IndirectIndexed), (Cpu6502::ex_adc, AddressMode::ZeroPageX), (Cpu6502::ex_adc, AddressMode::AbsoluteY), (Cpu6502::ex_adc, AddressMode::AbsoluteX) ],
    /*     4*/ [ (Cpu6502::ex_sta, AddressMode::IndexedIndirect), (Cpu6502::ex_sta, AddressMode::ZeroPage), (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_sta, AddressMode::Absolute), (Cpu6502::ex_sta, AddressMode::IndirectIndexed), (Cpu6502::ex_sta, AddressMode::ZeroPageX), (Cpu6502::ex_sta, AddressMode::AbsoluteY), (Cpu6502::ex_sta, AddressMode::AbsoluteX) ],
    /*     5*/ [ (Cpu6502::ex_lda, AddressMode::IndexedIndirect), (Cpu6502::ex_lda, AddressMode::ZeroPage), (Cpu6502::ex_lda, AddressMode::Immediate), (Cpu6502::ex_lda, AddressMode::Absolute), (Cpu6502::ex_lda, AddressMode::IndirectIndexed), (Cpu6502::ex_lda, AddressMode::ZeroPageX), (Cpu6502::ex_lda, AddressMode::AbsoluteY), (Cpu6502::ex_lda, AddressMode::AbsoluteX) ],
    /*     6*/ [ (Cpu6502::ex_cmp, AddressMode::IndexedIndirect), (Cpu6502::ex_cmp, AddressMode::ZeroPage), (Cpu6502::ex_cmp, AddressMode::Immediate), (Cpu6502::ex_cmp, AddressMode::Absolute), (Cpu6502::ex_cmp, AddressMode::IndirectIndexed), (Cpu6502::ex_cmp, AddressMode::ZeroPageX), (Cpu6502::ex_cmp, AddressMode::AbsoluteY), (Cpu6502::ex_cmp, AddressMode::AbsoluteX) ],
    /*     7*/ [ (Cpu6502::ex_sbc, AddressMode::IndexedIndirect), (Cpu6502::ex_sbc, AddressMode::ZeroPage), (Cpu6502::ex_sbc, AddressMode::Immediate), (Cpu6502::ex_sbc, AddressMode::Absolute), (Cpu6502::ex_sbc, AddressMode::IndirectIndexed), (Cpu6502::ex_sbc, AddressMode::ZeroPageX), (Cpu6502::ex_sbc, AddressMode::AbsoluteY), (Cpu6502::ex_sbc, AddressMode::AbsoluteX) ],
    ],

    // c = 2
    [          // b = 0                                        = 1                                       = 2                                          = 3                                       = 4                                       = 5                                        = 6                                       = 7
    /* a = 0*/ [ (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_asl, AddressMode::ZeroPage), (Cpu6502::ex_asl, AddressMode::Accumulator), (Cpu6502::ex_asl, AddressMode::Absolute), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_asl, AddressMode::ZeroPageX), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_asl, AddressMode::Implicit) ],
    /*     1*/ [ (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_rol, AddressMode::ZeroPage), (Cpu6502::ex_rol, AddressMode::Accumulator), (Cpu6502::ex_rol, AddressMode::Absolute), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_rol, AddressMode::ZeroPageX), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_rol, AddressMode::Implicit) ],
    /*     2*/ [ (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_lsr, AddressMode::ZeroPage), (Cpu6502::ex_lsr, AddressMode::Accumulator), (Cpu6502::ex_lsr, AddressMode::Absolute), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_lsr, AddressMode::ZeroPageX), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_lsr, AddressMode::Implicit) ],
    /*     3*/ [ (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_ror, AddressMode::ZeroPage), (Cpu6502::ex_ror, AddressMode::Accumulator), (Cpu6502::ex_ror, AddressMode::Absolute), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_ror, AddressMode::ZeroPageX), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_ror, AddressMode::Implicit) ],
    /*     4*/ [ (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_stx, AddressMode::ZeroPage), (Cpu6502::ex_txa, AddressMode::Implicit),    (Cpu6502::ex_stx, AddressMode::Absolute), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_stx, AddressMode::ZeroPageY), (Cpu6502::ex_txs, AddressMode::Implicit), (Cpu6502::op_err, AddressMode::Implicit) ],
    /*     5*/ [ (Cpu6502::ex_ldx, AddressMode::Immediate), (Cpu6502::ex_ldx, AddressMode::ZeroPage), (Cpu6502::ex_tax, AddressMode::Implicit),    (Cpu6502::ex_ldx, AddressMode::Absolute), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_ldx, AddressMode::ZeroPageY), (Cpu6502::ex_tsx, AddressMode::Implicit), (Cpu6502::ex_ldx, AddressMode::Implicit) ],
    /*     6*/ [ (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_dec, AddressMode::ZeroPage), (Cpu6502::ex_dex, AddressMode::Implicit),    (Cpu6502::ex_dec, AddressMode::Absolute), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_dec, AddressMode::ZeroPageX), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_dec, AddressMode::Implicit) ],
    /*     7*/ [ (Cpu6502::op_err, AddressMode::NA),        (Cpu6502::ex_inc, AddressMode::ZeroPage), (Cpu6502::ex_nop, AddressMode::Implicit),    (Cpu6502::ex_inc, AddressMode::Absolute), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_inc, AddressMode::ZeroPageX), (Cpu6502::op_err, AddressMode::NA),       (Cpu6502::ex_inc, AddressMode::Implicit) ],
    ],
];
