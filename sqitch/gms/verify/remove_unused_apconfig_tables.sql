-- Verify remove_unused_apconfig_tables
DO $$
BEGIN

IF EXISTS(SELECT * from information_schema.columns
    WHERE table_schema = 'config'
    AND table_name = 'analysis_project'
    AND column_name = 'analysis_menu_item_id') THEN
    RAISE EXCEPTION 'analyis_menu_item_id still exists!';
END IF;

IF EXISTS(SELECT * from information_schema.columns
    WHERE table_schema = 'config'
    AND table_name = 'analysis_project'
    AND column_name = 'configuration_set_id') THEN
    RAISE EXCEPTION 'configuration_set_id still exists!';
END IF;

IF EXISTS(SELECT * from information_schema.columns
    WHERE table_schema = 'config'
    AND table_name = 'analysis_menu_item') THEN
    RAISE EXCEPTION 'analysis_menu_item still exists!';
END IF;

IF EXISTS(SELECT * from information_schema.columns
    WHERE table_schema = 'config'
    AND table_name = 'set') THEN
    RAISE EXCEPTION 'config.set still exists!';
END IF;

END $$;
