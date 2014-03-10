-- Deploy add_bisulfite_column_to_fragment_library
-- requires: instrument_fragment_library

BEGIN;

DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM information_schema.columns
            WHERE table_schema='instrument'
                AND table_name='fragment_library'
                AND column_name='bisulfite_conversion'
        ) THEN
            ALTER TABLE instrument.fragment_library ADD COLUMN bisulfite_conversion VARCHAR(32);
        END IF;
    END;
$$;

COMMIT;
