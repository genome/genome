-- Verify drop_featurelist_file_id

DO $$
    BEGIN
        IF EXISTS(SELECT * FROM information_schema.columns
            WHERE table_schema = 'model'
            AND table_name = 'feature_list'
            AND column_name = 'file_id') THEN
            RAISE EXCEPTION 'file_id still exists!';
        END IF;
    END
$$
