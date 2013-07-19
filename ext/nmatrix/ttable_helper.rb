#!/usr/bin/ruby

# A helper file for generating and maintaining template tables.

def nullify(disabled) #:nodoc:
  DTYPES.map { |t| if disabled.include?(t) then :NULL else t end }
end

DTYPES = [
          :uint8_t,
          :int8_t,
          :int16_t,
          :int32_t,
          :int64_t,
          :float32_t,
          :float64_t,
          :'nm::Complex64',
          :'nm::Complex128',
          :'nm::Rational32',
          :'nm::Rational64',
          :'nm::Rational128',
          :'nm::RubyObject'
         ]

ITYPES = [
          :uint8_t,
          :uint16_t,
          :uint32_t,
          :uint64_t
         ]

EWOPS = [
         :'nm::EW_ADD',
         :'nm::EW_SUB',
         :'nm::EW_MUL',
         :'nm::EW_DIV',
         :'nm::EW_POW',
         :'nm::EW_MOD',
         :'nm::EW_EQEQ',
         :'nm::EW_NEQ',
         :'nm::EW_LT',
         :'nm::EW_GT',
         :'nm::EW_LEQ',
         :'nm::EW_GEQ'
        ]

LR_ALLOWED = {
  :uint8_t	 		=> nullify([:RubyObject]),
  :int8_t				=> nullify([:RubyObject]),
  :int16_t			=> nullify([:RubyObject]),
  :int32_t			=> nullify([:RubyObject]),
  :int64_t			=> nullify([:RubyObject]),
  :float32_t		=> nullify([:RubyObject]),
  :float64_t		=> nullify([:RubyObject]),
  :Complex64		=> nullify([:RubyObject]),
  :Complex128		=> nullify([:RubyObject]),
  :Rational32		=> nullify([:float32_t, :float64_t, :'nm::Complex64', :'nm::Complex128', :'nm::RubyObject']),
  :Rational64		=> nullify([:float32_t, :float64_t, :'nm::Complex64', :'nm::Complex128', :'nm::RubyObject']),
  :Rational128	=> nullify([:float32_t, :float64_t, :'nm::Complex64', :'nm::Complex128', :'nm::RubyObject']),
  :RubyObject		=> nullify(DTYPES - [:RubyObject])
}

lines =
  case ARGV[0]
  when 'OPLR'
    '{' +
      EWOPS.map do |op|

    '{' +
      DTYPES.map do |l_dtype|

      '{' +
        LR_ALLOWED[l_dtype].map do |r_dtype|
        if r_dtype == :NULL
          'NULL'
        else
          "fun<#{op}, #{l_dtype}, #{r_dtype}>"
        end
      end.join(', ') +
        '}'

    end.join(",\n") +
      '}'

  end.join(",\n") +
      '}'

  when 'OPID'
    '{' +
      EWOPS.map do |op|
    '{' +
      ITYPES.map do |itype|
      '{' +
        DTYPES.map do |dtype|

        if dtype == :NULL
          'NULL'
        else
          "fun<#{op}, #{itype}, #{dtype}>"
        end

      end.join(",") +
        '}'
    end.join(",\\\n") +
      '}'
  end.join(",\\\n") +
      '}'

  when 'LR'
    '{' +
      DTYPES.map do |l_dtype|

    '{' +
      LR_ALLOWED[l_dtype].map do |r_dtype|
      if r_dtype == :NULL
        'NULL'
      else
        "fun<#{l_dtype}, #{r_dtype}>"
      end
    end.join(', ') +
      '}'

  end.join(",\n") +
      '}'
  end

puts lines
