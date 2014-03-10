-- Deploy add_clia_fields_to_analysis_projects
-- requires: config_analysis_project

BEGIN;

DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM information_schema.columns
            WHERE table_schema='config'
                AND table_name='analysis_project'
                AND column_name='is_clia'
        ) THEN
            ALTER TABLE config.analysis_project ADD COLUMN is_clia BOOLEAN;
        END IF;
    END;
$$;

DO $$
    BEGIN
        IF NOT EXISTS (
            SELECT 1 FROM information_schema.columns
            WHERE table_schema='config'
                AND table_name='analysis_project'
                AND column_name='run_as'
        ) THEN
            ALTER TABLE config.analysis_project ADD COLUMN run_as VARCHAR(64);
        END IF;
    END;
$$;

COMMIT;
