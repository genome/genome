-- Deploy model.model.created_by
-- requires: model_model

BEGIN;

DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM information_schema.columns
            WHERE table_schema='model'
                AND table_name='model'
                AND column_name='created_by'
        ) THEN
            ALTER TABLE model.model ADD COLUMN created_by VARCHAR(64);
        END IF;
    END;
$$;

COMMIT;
